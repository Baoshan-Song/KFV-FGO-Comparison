import asyncio
from copy import deepcopy
from pathlib import Path
from tempfile import TemporaryDirectory
from threading import get_ident
from time import sleep, monotonic
import unittest
from unittest.mock import patch

import numpy as np
from colab_ui import ComparisonUI
from colab_comparison import prepare_comparison
from live_comparison import LiveComparison
from run_cache import load_result_cache, method_key
from notebook_transport import ColabPanel, close_panels


class ColabRequestTests(unittest.TestCase):
    """No event loop, no await and no cell kept running to drive the UI."""
    def setUp(self):
        self.callbacks={}
        class Output:
            register_callback=lambda _,name,fn:self.callbacks.__setitem__(name,fn)
            unregister_callback=lambda _,name:self.callbacks.pop(name)
        self.output=Output()

    def tearDown(self):
        close_panels()

    def panel(self,module='KFV vs FGO'):
        ui=ComparisonUI(module);ui.sim['num_steps'].value=30
        panel=ColabPanel(ui,self.output)
        with patch('notebook_transport.display'):
            panel._ipython_display_()
        self.assertTrue(ui.run.disabled)
        self.callbacks[panel.name]()
        self.assertFalse(ui.run.disabled)
        self.assertIn('Colab connected',ui.connection.value)
        return ui,panel

    def finish(self,ui,panel):
        deadline=monotonic()+5
        while ui._busy and monotonic()<deadline:
            self.callbacks[panel.name]()
            sleep(.01)
        self.assertFalse(ui._busy,ui.status.value)
        self.assertFalse(ui.run.disabled)
        self.assertTrue(ui.stop.disabled)

    def test_colab_runs_and_restarts_without_notebook_asyncio(self):
        ui,panel=self.panel()
        for _ in range(2):
            with patch('colab_ui.asyncio.get_running_loop',side_effect=AssertionError('Colab must not use asyncio')):
                ui.run.click()
            self.finish(ui,panel)
            self.assertIn('Complete',ui.status.value)
            self.assertFalse(ui.details.disabled)
            self.assertFalse(hasattr(ui,'ui_task'))
            state=panel.poll()
            self.assertIn('class="play"',state['viewer']['html'])
            self.assertIn('requestAnimationFrame',state['viewer']['html'])
            self.assertEqual(ui.live_container.layout.display,'none')
            self.assertNotIn('viewer',panel.poll(viewer_revision=state['viewer']['revision']))

    def test_stop_is_acknowledged_then_can_run_again(self):
        ui,panel=self.panel()
        def loading(*args,**kwargs):
            while not kwargs['stop'].is_set():sleep(.005)
            kwargs['progress']('Reading observations')
        with patch('colab_ui.prepare_comparison',side_effect=loading):
            ui.run.click()
            ui.stop.click()
            self.assertIn('Stopping',ui.status.value)
            self.finish(ui,panel)
        self.assertIn('Stopped',ui.status.value)
        ui.run.click();self.finish(ui,panel)
        self.assertIn('Complete',ui.status.value)

    def test_cached_real_result_uses_same_request_path(self):
        ui,panel=self.panel()
        ui.source.value='Real GNSS/INS'
        with patch('colab_ui.prepare_comparison', wraps=prepare_comparison) as prepare:
            ui.run.click();self.finish(ui,panel)
        self.assertIn('Loaded saved results',ui.status.value)
        self.assertTrue(all(call.kwargs.get('load_observations') is False for call in prepare.call_args_list))
        self.assertEqual(ui._live.timeline.max,656)
        self.assertIn('class="play"',ui._viewer_html)

    def test_both_sliding_window_defaults_have_saved_results(self):
        ui,panel=self.panel('FGO vs FGO')
        ui.source.value='Real GNSS/INS'
        ui.run.click();self.finish(ui,panel)
        self.assertIn('Loaded saved results',ui.status.value)
        self.assertEqual(ui.cached_payload['names'],['SW-FGO (W=5)','SW-FGO (W=20)'])

    def test_startup_failure_visible_and_retry_enabled(self):
        ui,panel=self.panel()
        with patch.object(ui.left,'values',side_effect=ValueError('bad configuration')):
            ui.run.click()
        self.assertIn('Start failed: ValueError: bad configuration',ui.status.value)
        self.assertFalse(ui.run.disabled)
        self.assertFalse(ui._busy)
        ui.run.click();self.finish(ui,panel)
        self.assertIn('Complete',ui.status.value)

    def test_worker_failure_visible(self):
        ui,panel=self.panel()
        with patch('colab_ui.prepare_comparison',side_effect=ValueError('missing observations')):
            ui.run.click();self.finish(ui,panel)
        self.assertIn('missing observations',ui.status.value)

    def test_reopening_module_closes_old_callback_only(self):
        ui,old=self.panel()
        other_ui,other=self.panel('KFV vs KFV')
        _,new=self.panel()
        self.assertTrue(old.closed)
        self.assertTrue(ui.run.disabled)
        self.assertNotIn(old.name,self.callbacks)
        self.assertIn(new.name,self.callbacks)
        self.assertIn(other.name,self.callbacks)
        self.assertFalse(other_ui.run.disabled)


class NotebookJobTests(unittest.IsolatedAsyncioTestCase):
    async def test_loading_remains_responsive_and_widgets_stay_on_main_thread(self):
        ui=ComparisonUI('KFV vs FGO');ui.sim['num_steps'].value=30
        main=get_ident();threads=[];ticks=[]
        class TrackedView(LiveComparison):
            def __init__(self,*args):
                threads.append(get_ident());super().__init__(*args)
            def receive(self,*args):
                threads.append(get_ident());super().receive(*args)
        def delayed(*args,**kwargs):
            self.assertNotEqual(get_ident(),main)
            kwargs['progress']('Reading test observations')
            sleep(.15)
            return prepare_comparison(*args,**kwargs)
        with patch('colab_ui.prepare_comparison',side_effect=delayed),patch('colab_ui.LiveComparison',TrackedView):
            ui._run(None)
            while not ui.ui_task.done():
                ticks.append(ui.status.value);await asyncio.sleep(.02)
            await ui.ui_task
        self.assertTrue(any('Reading test observations' in t for t in ticks))
        self.assertTrue(threads and all(t==main for t in threads))
        self.assertIn('Complete',ui.status.value)
        self.assertFalse(ui.run.disabled)
        self.assertFalse(ui.details.disabled)

    async def test_stop_during_loading_restores_controls(self):
        ui=ComparisonUI('KFV vs KFV')
        def loading(*args,**kwargs):
            while not kwargs['stop'].is_set():sleep(.005)
            kwargs['progress']('Loading')
        with patch('colab_ui.prepare_comparison',side_effect=loading):
            ui._run(None);await asyncio.sleep(.03);ui.stop_event.set()
            await asyncio.wait_for(ui.ui_task,2)
        self.assertIn('Stopped',ui.status.value)
        self.assertFalse(ui.run.disabled)
        self.assertTrue(ui.stop.disabled)

    async def test_loading_failure_is_visible_not_stuck(self):
        ui=ComparisonUI('KFV vs FGO')
        with patch('colab_ui.prepare_comparison',side_effect=ValueError('test missing observations')):
            ui._run(None);await asyncio.wait_for(ui.ui_task,2)
        self.assertIn('test missing observations',ui.status.value)
        self.assertFalse(ui.run.disabled)


class CacheTests(unittest.TestCase):
    def test_configuration_change_cannot_reuse_saved_run(self):
        _,a,b=prepare_comparison('KFV vs FGO','Real GNSS/INS',
            dict(method='iEKF'),dict(method='Re-FGO'),load_observations=False)
        changed=deepcopy(a);changed['config'].kfv.q[0,0]*=2
        self.assertNotEqual(method_key(a,'test-input'),method_key(changed,'test-input'))
        self.assertNotEqual(method_key(a,'test-input'),method_key(a,'changed-input'))
        self.assertNotEqual(method_key(a,'test-input'),method_key(b,'test-input'))

    def test_saved_replay_keeps_signed_up_and_missing_residuals(self):
        n=4;trajectory=[[0,1,2,3],[0,0,0,0],[0,1,2,3]]
        payload=dict(real=True,names=['A','B'],kinds=['kfv','fgo'],frames=n,dt=1.,time=[0,1,2,4],
            trajectories=[trajectory,trajectory],truth=trajectory,has_truth=True,
            errors=[dict(horizontal=[None,2,3,4],up=[None,-1,-2,-3],**{'3d':[None,3,4,5]})]*2,
            residuals={mode:[[[None,2,None,3],[None,None,4,None]]]*2 for mode in ('postfit','prefit')})
        view=LiveComparison.from_payload(payload)
        self.assertEqual(view.timeline.max,3)
        self.assertEqual(view.time[-1],4.)
        np.testing.assert_allclose(view.norms[0,0,1:],[2,4,3])
        view.view.value='Position error';view.component.value='Up'
        np.testing.assert_allclose(view.scalar_errors['up'][0,1:],[-1,-2,-3])

    def test_corrupt_cache_falls_back_to_compute(self):
        _,a,b=prepare_comparison('KFV vs FGO','Real GNSS/INS',
            dict(method='iEKF'),dict(method='Re-FGO'),load_observations=False)
        with TemporaryDirectory() as directory,patch('run_cache.input_fingerprint',return_value='input'):
            (Path(directory)/(method_key(a,'input')+'.json.gz')).write_bytes(b'incomplete gzip')
            self.assertIsNone(load_result_cache([a,b],[directory]))


if __name__=='__main__':unittest.main()
