import unittest

from colab_ui import ComparisonUI


class TestSettingsSections(unittest.TestCase):
    def test_simulation_settings_open_and_preserve_shared_controls(self):
        ui = ComparisonUI('KFV vs KFV')
        groups = ui.sim_panel.children
        self.assertEqual(len(groups), 3)
        for group in groups:
            button, panel = group.children
            self.assertEqual(panel.layout.display, 'none')
            button.click()
            self.assertEqual(panel.layout.display, '')
            button.click()
            self.assertEqual(panel.layout.display, 'none')
        # Grouping must retain the same model objects consumed by the run path.
        self.assertIs(groups[0].children[1], ui.sim_data)
        ui.sim['num_steps'].value = 150
        ui.sim['outlier_weight'].value = .35
        ui.sim['initial_error_x'].value = 25
        self.assertEqual(ui.sim['num_steps'].value, 150)
        self.assertAlmostEqual(ui.sim['outlier_weight'].value, .35)
        self.assertEqual(ui.sim['initial_error_x'].value, 25)

    def test_source_switch_keeps_both_groups_available(self):
        for module in ('KFV vs KFV', 'KFV vs FGO', 'FGO vs FGO'):
            ui = ComparisonUI(module)
            ui.source.value = 'Real GNSS/INS'
            self.assertEqual(ui.sim_panel.layout.display, 'none')
            self.assertEqual(ui.real_panel.layout.display, '')
            button, panel = ui.real_panel.children[0].children
            button.click()
            self.assertEqual(panel.layout.display, '')
            ui.source.value = 'Simulation'
            self.assertEqual(ui.sim_panel.layout.display, '')
            self.assertEqual(ui.real_panel.layout.display, 'none')


if __name__ == '__main__':
    unittest.main()
