"""Explicit Colab requests for live updates after a widget callback returns.

Only numerical work runs on the worker. Widgets are created and updated inside
short kernel requests from Colab's documented register_callback/invokeFunction
API, so progress does not depend on a background notebook asyncio task.
"""
from html import escape
import json
import base64
import gzip
from pathlib import Path
from uuid import uuid4

from IPython.display import HTML, JSON, display

UI_REVISION = '2026-11-22.1'
_PANELS = {}


def encode_viewer(payload):
    text=json.dumps(payload,allow_nan=False,separators=(',',':'))
    if len(text)<64000:return dict(encoding='json',data=payload)
    return dict(encoding='gzip-base64',data=base64.b64encode(gzip.compress(text.encode(),mtime=0)).decode())


def colab_output():
    try:
        from google.colab import output
    except ImportError:
        return None
    return output


def close_panels():
    for panel in list(_PANELS.values()):
        panel.close()


class ColabPanel:
    def __init__(self, ui, output):
        self.ui = ui
        self.output = output
        self.name = 'kfv_progress_' + uuid4().hex
        self.connected = False
        self.closed = False
        self.mounted = False
        ui._colab_bridge = self
        ui.run.disabled = True
        ui.connection.value = f'<small>UI {UI_REVISION} · Connecting Colab controls…</small>'

    def poll(self, action='poll', viewer_revision=-1, rendered_revision=-1):
        if self.closed:
            return dict(closed=True, active=False)
        if not self.connected:
            self.connected = True
            self.ui.run.disabled = False
            self.ui.connection.value = f'<small>UI {UI_REVISION} · Colab connected</small>'
        if action == 'stop':
            self.ui._stop()
        self.ui._viewer_ready(rendered_revision)
        self.ui._drain_events()
        state=dict(closed=False, active=self.ui._busy)
        if viewer_revision!=self.ui._viewer_revision:
            state['viewer']=dict(revision=self.ui._viewer_revision,packet=self.ui._viewer_packet)
        return state

    def _callback(self, action='poll', viewer_revision=-1, rendered_revision=-1):
        return JSON(self.poll(action,viewer_revision,rendered_revision))

    def close(self):
        self.closed = True
        self.ui._closed = True
        self.ui._stop()
        self.ui.run.disabled = True
        self.ui.connection.value = f'<small>UI {UI_REVISION} · Replaced by a new panel</small>'
        if self.mounted:
            self.output.unregister_callback(self.name)
        if _PANELS.get(self.ui.module) is self:
            del _PANELS[self.ui.module]

    def _ipython_display_(self):
        if not self.mounted:
            previous = _PANELS.get(self.ui.module)
            if previous is not None and previous is not self:
                previous.close()
            _PANELS[self.ui.module] = self
            self.output.register_callback(self.name, self._callback)
            self.mounted = True
        display(self.ui.widget)
        display(HTML(self._script()))

    def _script(self):
        # Mount the original renderer in this Colab output, not in a second
        # sandboxed iframe. It receives data without injecting executable HTML.
        template=(Path(__file__).parent/'kfv_fgo'/'comparison.html').read_text()
        plot=template.replace('__ID__',self.name+'_plot').replace('__PAYLOAD__','null').replace('__INTERVAL__','80')
        bridge=r'''<div id="__ID__" style="font:13px Arial;color:#52677b" role="status"></div>
<script>
(() => {
  const note=document.getElementById('__ID__');
  const plot=document.getElementById('__ID___plot');
  const name=__NAME__;
  let closed=false, viewerRevision=-1, renderedRevision=-1, plotError='', failures=0;
  plot.addEventListener('comparison-ready',event=>{
    renderedRevision=event.detail.revision;
    plotError='';
    try { google.colab.output?.setIframeHeight?.(-1,true,{interactive:true,maxHeight:99999}); } catch (_) {}
  });
  plot.addEventListener('comparison-error',event=>{
    plotError='Cannot display the full plot: '+event.detail.message+'. Try Run comparison again.';
    note.textContent=plotError;
  });
  async function showViewer(viewer) {
    let data=null;
    try {
      if(viewer.packet?.encoding==='json')data=viewer.packet.data;
      else if(viewer.packet?.encoding==='gzip-base64'){
        const bytes=Uint8Array.from(atob(viewer.packet.data),c=>c.charCodeAt(0));
        data=JSON.parse(await new Response(new Blob([bytes]).stream().pipeThrough(new DecompressionStream('gzip'))).text());
      }else if(viewer.packet)throw new Error('Unknown plot data format');
      plot.dispatchEvent(new CustomEvent('comparison-data',{detail:{revision:viewer.revision,data}}));
      if(renderedRevision!==viewer.revision)throw new Error('The renderer did not confirm the result');
      viewerRevision=viewer.revision;
    }catch(error){plotError='Cannot load the full plot: '+String(error)+'. Try Run comparison again.';}
  }
  async function poll() {
    if (closed || !note.isConnected) return;
    const warning=setTimeout(() => {
      note.textContent='Colab is not responding. Check the runtime connection.';
    }, 10000);
    try {
      const response=await google.colab.kernel.invokeFunction(name, [], {
        viewer_revision:viewerRevision,rendered_revision:renderedRevision});
      const state=response.data['application/json'];
      clearTimeout(warning);
      if (!state || typeof state.active !== 'boolean') throw new Error('Invalid progress response');
      failures=0;
      closed=state.closed;
      if(state.viewer)await showViewer(state.viewer);
      note.textContent=plotError;
      if (!closed) setTimeout(poll, state.active ? 200 : 1000);
    } catch (error) {
      clearTimeout(warning);
      failures+=1;
      note.textContent='Reconnecting to Colab… If this continues, check the runtime connection.';
      if(!closed && note.isConnected)setTimeout(poll,Math.min(5000,500*failures));
    }
  }
  poll();
})();
</script>'''
        return plot+bridge.replace('__ID__',escape(self.name)).replace('__NAME__',json.dumps(self.name))
