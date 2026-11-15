"""Explicit Colab requests for live updates after a widget callback returns.

Only numerical work runs on the worker. Widgets are created and updated inside
short kernel requests from Colab's documented register_callback/invokeFunction
API, so progress does not depend on a background notebook asyncio task.
"""
from html import escape
import json
from uuid import uuid4

from IPython.display import HTML, JSON, display

UI_REVISION = '2026-11-15.2'
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

    def poll(self, action='poll'):
        if self.closed:
            return dict(closed=True, active=False)
        if not self.connected:
            self.connected = True
            self.ui.run.disabled = False
            self.ui.connection.value = f'<small>UI {UI_REVISION} · Colab connected</small>'
        if action == 'stop':
            self.ui._stop()
        self.ui._drain_events()
        return dict(closed=False, active=self.ui._busy)

    def _callback(self, action='poll'):
        return JSON(self.poll(action))

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
        return '''<div id="__ID__" style="font:13px Arial;color:#52677b" role="status"></div>
<script>
(() => {
  const note=document.getElementById('__ID__');
  const name=__NAME__;
  let closed=false;
  async function poll() {
    if (closed) return;
    const warning=setTimeout(() => {
      note.textContent='Colab is not responding. Check the runtime connection; progress resumes after it reconnects.';
    }, 10000);
    try {
      const response=await google.colab.kernel.invokeFunction(name, [], {});
      const state=response.data['application/json'];
      clearTimeout(warning);
      if (!state || typeof state.active !== 'boolean') throw new Error('Invalid progress response');
      closed=state.closed;
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
