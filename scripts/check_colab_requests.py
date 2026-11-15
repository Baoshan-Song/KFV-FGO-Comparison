"""Exercise notebook buttons and Colab-style polling in separate kernel requests.

No cell awaits the comparison task. The kernel returns idle after every request.
The browser transport is emulated; this does not claim a Google-hosted UI test.
Requires jupyter_client and ipykernel (development only).
"""
from pathlib import Path
from time import monotonic, sleep

from jupyter_client import KernelManager


def main():
    manager=KernelManager(kernel_name='python3')
    manager.start_kernel(cwd=str(Path(__file__).resolve().parents[1]))
    client=manager.blocking_client();client.start_channels();client.wait_for_ready(timeout=30)
    def execute(code):
        ident=client.execute(code,store_history=False)
        result=None
        while True:
            message=client.get_iopub_msg(timeout=30)
            if message['parent_header'].get('msg_id')!=ident:continue
            content=message['content']
            if message['msg_type']=='error':raise RuntimeError(content)
            if message['msg_type'] in ('execute_result','display_data'):
                result=content['data'].get('application/json',result)
            if message['msg_type']=='status' and content['execution_state']=='idle':return result

    def click(comm_id):
        client.session.send(client.shell_channel.socket,'comm_msg',dict(
            comm_id=comm_id,data=dict(method='custom',content=dict(event='click'))))

    def poll_until_done():
        deadline=monotonic()+30
        while monotonic()<deadline:
            state=execute("panel._callback()")
            if not state['active']:return
            sleep(.1)
        raise AssertionError('Comparison did not complete')

    try:
        execute('''from colab_ui import ComparisonUI
from notebook_transport import ColabPanel
from IPython.display import display,JSON
from unittest.mock import patch
callbacks={}
class Output:
    def register_callback(self,name,fn):callbacks[name]=fn
    def unregister_callback(self,name):del callbacks[name]
''')
        for source,reuse in [('Simulation',True),('Real GNSS/INS',True),('Real GNSS/INS',False)]:
            ids=execute(f'''ui=ComparisonUI('KFV vs FGO')
ui.sim['num_steps'].value=300
ui.source.value={source!r}
ui.use_cache.value={reuse!r}
panel=ColabPanel(ui,Output())
display(panel)
panel.poll()
JSON(dict(run=ui.run.model_id,stop=ui.stop.model_id))''')
            click(ids['run']);poll_until_done()
            state=execute("JSON(dict(status=ui.status.value,enabled=not ui.run.disabled,async_task=hasattr(ui,'ui_task'),timeline=ui._live.timeline.max))")
            assert state['enabled'] and not state['async_task'],state
            assert ('Loaded saved results' if source!='Simulation' and reuse else 'Complete') in state['status'],state
            assert state['timeline']==(299 if source=='Simulation' else 656),state
            print('PASS',source,'cache' if reuse else 'fresh',state,flush=True)

        ids=execute('''import time
from colab_comparison import prepare_comparison
ui=ComparisonUI('KFV vs KFV')
panel=ColabPanel(ui,Output());display(panel);panel.poll()
def slow_loading(*args,**kwargs):
    while not kwargs['stop'].is_set():time.sleep(.01)
    kwargs['progress']('Loading')
loading_patch=patch('colab_ui.prepare_comparison',side_effect=slow_loading)
loading_patch.start()
JSON(dict(run=ui.run.model_id,stop=ui.stop.model_id))''')
        click(ids['run'])
        assert execute('panel._callback()')['active']
        click(ids['stop']);poll_until_done()
        state=execute("loading_patch.stop();JSON(dict(status=ui.status.value,enabled=not ui.run.disabled))")
        assert state['enabled'] and 'Stopped' in state['status'],state
        click(ids['run']);poll_until_done()
        state=execute("JSON(dict(status=ui.status.value,enabled=not ui.run.disabled))")
        assert state['enabled'] and 'Complete' in state['status'],state
        print('PASS actual Stop button, then Run again',state,flush=True)
    finally:
        client.stop_channels();manager.shutdown_kernel(now=True)


if __name__=='__main__':main()
