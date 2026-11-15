import unittest
from unittest.mock import patch
import numpy as np

from colab_ui import ComparisonUI
from colab_comparison import prepare_comparison


class TestSettingsSections(unittest.TestCase):
    def test_primary_settings_stay_visible_and_extras_expand(self):
        ui = ComparisonUI('KFV vs KFV')
        self.assertNotEqual(ui.sim_data.layout.display, 'none')
        self.assertNotEqual(ui.sim_faults.layout.display, 'none')
        button, panel = ui.sim_extra.children[0].children
        self.assertIs(panel, ui.sim_initial)
        self.assertEqual(panel.layout.display, 'none')
        button.click()
        self.assertEqual(panel.layout.display, '')
        button.click()
        self.assertEqual(panel.layout.display, 'none')
        self.assertIn('kfv-card-data', ui.data_card._dom_classes)

    def test_source_switch_preserves_edited_values(self):
        for module in ('KFV vs KFV', 'KFV vs FGO', 'FGO vs FGO'):
            ui = ComparisonUI(module)
            ui.sim['anchor_radius'].value = 105
            ui.sim['outlier_weight'].value = .35
            ui.sim['outlier_mean'].value = 6
            ui.sim['outlier_sigma'].value = 15
            ui.source.value = 'Real GNSS/INS'
            ui.real['minimum_snr_dbhz'].value = 25
            self.assertEqual(ui.sim_panel.layout.display, 'none')
            self.assertEqual(ui.real_panel.layout.display, '')
            button, panel = ui.real_extra.children[0].children
            button.click()
            self.assertEqual(panel.layout.display, '')
            ui.source.value = 'Simulation'
            self.assertEqual(ui.sim_panel.layout.display, '')
            self.assertEqual(ui.real_panel.layout.display, 'none')
            self.assertEqual(ui.sim['anchor_radius'].value, 105)
            self.assertEqual(ui.sim['outlier_weight'].value, .35)
            self.assertEqual(ui.sim['outlier_mean'].value, 6)
            self.assertEqual(ui.sim['outlier_sigma'].value, 15)
            ui.source.value = 'Real GNSS/INS'
            self.assertEqual(ui.real['minimum_snr_dbhz'].value, 25)

    def test_edited_simulation_values_reach_run_request(self):
        ui = ComparisonUI('KFV vs KFV')
        ui._colab_bridge = object()
        values = dict(anchor_radius=105.,outlier_weight=.4,outlier_mean=5.,
                      outlier_sigma=12.,white_sigma=.3,seed=11,num_steps=60)
        for key, value in values.items():ui.sim[key].value = value
        with patch.object(ui, '_start_worker') as worker:ui.run.click()
        self.assertTrue(worker.called)
        options = worker.call_args.args[1]
        for key, value in values.items():self.assertEqual(options['simulation'][key], value)
        self.assertEqual(ui.progress.layout.display, '')
        ui._restore_controls()
        self.assertEqual(ui.progress.layout.display, 'none')

    def test_real_controls_configure_both_methods(self):
        ui = ComparisonUI('KFV vs FGO')
        ui.source.value = 'Real GNSS/INS'
        values = dict(minimum_elevation_deg=10.,minimum_snr_dbhz=25.,
                      position_process_std_m=.4,velocity_process_std_mps=.2,
                      bias_process_std_mps2=.03,clock_process_std_m=80.)
        for key,value in values.items():ui.real[key].value=value
        _,a,b=prepare_comparison(ui.module,ui.source.value,ui.left.values(),ui.right.values(),
            real={k:v.value for k,v in ui.real.items()},match_filter_settings=True,load_observations=False)
        for method in (a,b):
            cfg=method['config']
            self.assertEqual(cfg.gnss['minimum_elevation_deg'],10.)
            self.assertEqual(cfg.gnss['minimum_snr'],25.)
            expected=np.array([.4]*3+[.2]*3+[.03]*3+[80.])**2
            np.testing.assert_allclose(np.diag(getattr(cfg,method['kind']).q),expected)


if __name__ == '__main__':
    unittest.main()
