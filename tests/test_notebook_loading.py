import asyncio
from copy import deepcopy
from pathlib import Path
from tempfile import TemporaryDirectory
from threading import get_ident
from time import sleep
import unittest
from unittest.mock import patch

import numpy as np
from colab_ui import ComparisonUI
from colab_comparison import prepare_comparison
from live_comparison import LiveComparison
from run_cache import load_result_cache, method_key


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
