import unittest
from unittest.mock import patch
import numpy as np

from colab_ui import ComparisonUI
from colab_comparison import prepare_comparison


class ColabModulesTests(unittest.TestCase):
    def test_filter_controls_follow_effective_mode(self):
        ui=ComparisonUI('KFV vs KFV')
        c=ui.left
        self.assertEqual((c.iterations.value,c.kernel.value),(1,'none'))
        self.assertTrue(c.iterations.disabled and c.kernel.disabled and c.delta.disabled)
        c.method.value='iEKF'
        self.assertFalse(c.iterations.disabled)
        self.assertTrue(c.kernel.disabled)
        c.iterations.value=8
        c.method.value='rEKF'
        self.assertEqual(c.iterations.value,1)
        self.assertFalse(c.kernel.disabled)
        c.method.value='riEKF'
        self.assertEqual(c.iterations.value,8)
        self.assertFalse(c.iterations.disabled or c.kernel.disabled)

    def test_linked_graph_and_source_visibility(self):
        ui=ComparisonUI('KFV vs FGO')
        ui.left.method.value='EKF'
        self.assertEqual((ui.right.iterations.value,ui.right.kernel.value),(1,'none'))
        self.assertTrue(ui.right.iterations.disabled and ui.right.kernel.disabled)
        ui.source.value='Real GNSS/INS'
        self.assertEqual(ui.sim_panel.layout.display,'none')
        self.assertNotEqual(ui.real_panel.layout.display,'none')
        self.assertEqual(ui.right.numerical_row.layout.display,'none')
        ui.right.method.value='SW-FGO'
        self.assertFalse(ui.match.value)
        self.assertTrue(ui.match.disabled)
        self.assertFalse(ui.right.window.disabled)
        self.assertFalse(ui.right.iterations.disabled)
        ui.right.method.value='Re-FGO'
        self.assertEqual(ui.right.window.value,1)
        self.assertTrue(ui.right.window.disabled)

    def test_all_modules_share_inputs_and_noise_in_both_sources(self):
        for module,methods in [('KFV vs KFV',('EKF','riEKF')),
                               ('KFV vs FGO',('iEKF','Re-FGO')),
                               ('FGO vs FGO',('Re-FGO','SW-FGO'))]:
            for source in ('Simulation','Real GNSS/INS'):
                with self.subTest(module=module,source=source), patch('kfv_fgo.data.gnss_ins.GnssImuDataset'):
                    data,a,b=prepare_comparison(module,source,
                        dict(method=methods[0],iterations=3,kernel='none'),
                        dict(method=methods[1],iterations=4,kernel='huber',window=7),
                        simulation=dict(num_steps=8,white_sigma=.7),
                        real=dict(position_process_std_m=2.,velocity_process_std_mps=.4),
                        match_filter_settings=module=='KFV vs FGO')
                    x=getattr(a['config'],a['kind']);y=getattr(b['config'],b['kind'])
                    np.testing.assert_array_equal(x.q,y.q)
                    np.testing.assert_array_equal(x.p0,y.p0)
                    if source=='Simulation':self.assertAlmostEqual(x.r,.49)
                    else:self.assertAlmostEqual(x.q[0,0],4.)
                    if module=='KFV vs FGO':
                        self.assertEqual(y.max_iteration,3)
                        self.assertEqual(y.robust_kernel,'none')
                        self.assertTrue(y.imitate_kfv)
                    if module=='FGO vs FGO':self.assertEqual(y.window_size,7)


if __name__=='__main__':unittest.main()
