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
        self.assertFalse(ui.match_filter_settings)
        self.assertFalse(ui.right.window.disabled)
        self.assertFalse(ui.right.iterations.disabled)
        ui.right.method.value='Re-FGO'
        self.assertEqual(ui.right.window.value,1)
        self.assertTrue(ui.right.window.disabled)
        self.assertTrue(ui.match_filter_settings)
        self.assertTrue(ui.right.iterations.disabled and ui.right.kernel.disabled)

    def test_kernel_scales_are_remembered_and_matched(self):
        ui=ComparisonUI('KFV vs FGO')
        ui.left.method.value='riEKF'
        ui.left.kernel.value='tukey'
        self.assertAlmostEqual(ui.left.delta.value,4.685)
        self.assertEqual(ui.left.delta.value,ui.right.delta.value)
        ui.left.delta.value=6.
        ui.left.kernel.value='cauchy'
        self.assertAlmostEqual(ui.left.delta.value,2.385)
        ui.left.kernel.value='tukey'
        self.assertEqual(ui.left.delta.value,6.)
        self.assertEqual(ui.right.delta.value,6.)

    def test_recursive_graph_always_tracks_the_effective_filter(self):
        ui=ComparisonUI('KFV vs FGO')
        self.assertFalse(hasattr(ui,'match'))
        for mode in ('EKF','rEKF','iEKF','riEKF'):
            ui.left.method.value=mode
            if not ui.left.iterations.disabled:ui.left.iterations.value=7
            if not ui.left.kernel.disabled:ui.left.kernel.value='cauchy'
            for field in ('iterations','kernel','delta','continuation'):
                left,right=getattr(ui.left,field),getattr(ui.right,field)
                self.assertEqual(left.value,right.value,(mode,field))
                self.assertTrue(right.disabled,(mode,field))
            self.assertTrue(ui.right.numerical.disabled)
            self.assertEqual(ui.right.window.value,1)
            ui.right.method.value='SW-FGO'
            self.assertFalse(ui.right.iterations.disabled)
            ui.right.method.value='Re-FGO'
            self.assertEqual(ui.right.iterations.value,ui.left.iterations.value)
            self.assertTrue(ui.right.iterations.disabled)

    def test_graph_module_only_offers_independent_sliding_windows(self):
        ui=ComparisonUI('FGO vs FGO')
        self.assertFalse(ui.match_filter_settings)
        self.assertEqual((ui.left.window.value,ui.right.window.value),(5,20))
        for card in (ui.left,ui.right):
            self.assertEqual(tuple(card.method.options),('SW-FGO',))
            self.assertTrue(card.method.disabled)
            self.assertFalse(card.window.disabled)
        ui.left.window.value=10
        self.assertEqual(ui.right.window.value,20)

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
