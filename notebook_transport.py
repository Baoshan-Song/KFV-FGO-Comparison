"""Explicit Colab requests for live updates after a widget callback returns.

Only numerical work runs on the worker. Widgets are created and updated inside
short kernel requests from Colab's documented register_callback/invokeFunction
API, so progress does not depend on a background notebook asyncio task.
"""
from html import escape
import json
from uuid import uuid4

from IPython.display import HTML, JSON, display

UI_REVISION = '2026-11-15.3'
_PANELS = {}


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

    def poll(self, action='poll', viewer_revision=-1):
        if self.closed:
            return dict(closed=True, active=False)
        if not self.connected:
            self.connected = True
            self.ui.run.disabled = False
            self.ui.connection.value = f'<small>UI {UI_REVISION} · Colab connected</small>'
        if action == 'stop':
            self.ui._stop()
        self.ui._drain_events()
        state=dict(closed=False, active=self.ui._busy)
        if viewer_revision!=self.ui._viewer_revision:
            state['viewer']=dict(revision=self.ui._viewer_revision,html=self.ui._viewer_html)
        return state

    def _callback(self, action='poll', viewer_revision=-1):
        return JSON(self.poll(action,viewer_revision))

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
        # This is an ordinary Colab HTML output, not widgets.HTML (which strips
        # scripts). One request is in flight at a time, including when delayed.
        return r'''<div id="__ID__" style="font:13px Arial;color:#52677b" role="status"></div>
<iframe id="__ID___replay" title="Comparison replay" hidden
  sandbox="allow-scripts allow-downloads" style="width:100%;height:760px;border:0"></iframe>
<script>
(() => {
  const note=document.getElementById('__ID__');
  const replay=document.getElementById('__ID___replay');
  const name=__NAME__;
  let closed=false, viewerRevision=-1;
  window.addEventListener('message',event=>{
    if(event.source===replay.contentWindow && event.data?.channel===name &&
       Number.isFinite(event.data.height)) replay.style.height=Math.max(200,event.data.height)+'px';
  });
  function showReplay(viewer) {
    viewerRevision=viewer.revision;
    replay.hidden=!viewer.html;
    if(!viewer.html){replay.srcdoc='';return;}
    const resize=`<script>new ResizeObserver(()=>parent.postMessage({channel:${JSON.stringify(name)},height:Math.ceil(document.body.getBoundingClientRect().height)+8},'*')).observe(document.body);<\/script>`;
    replay.srcdoc='<!doctype html><meta charset="utf-8"><style>body{margin:0}</style>'+viewer.html+resize;
  }
  async function poll() {
    if (closed) return;
    const warning=setTimeout(() => {
      note.textContent='Colab is not responding. Check the runtime connection; progress resumes after it reconnects.';
    }, 10000);
    try {
      const response=await google.colab.kernel.invokeFunction(name, [], {viewer_revision:viewerRevision});
      const state=response.data['application/json'];
      clearTimeout(warning);
      if (!state || typeof state.active !== 'boolean') throw new Error('Invalid progress response');
      closed=state.closed;
      if(state.viewer)showReplay(state.viewer);
      note.textContent='';
      if (!closed) setTimeout(poll, state.active ? 200 : 1000);
    } catch (error) {
      clearTimeout(warning);
      closed=true;
      note.textContent='Live connection interrupted. Run this comparison cell again. '+String(error);
    }
  }
  poll();
})();
</script>'''.replace('__ID__', escape(self.name)).replace('__NAME__', json.dumps(self.name))
