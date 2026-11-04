import json
from pathlib import Path
import unittest
import numpy as np
from experiment_common import (make_data, make_config, run_estimator,
    measurement_residuals, metrics, comparison_payload, animate_pair)
from config.config import range_measurement


class ResidualTests(unittest.TestCase):
    def test_online_residual_is_not_smoothed_or_truth_error(self):
        data = make_data(num_steps=12)
        result, elapsed = run_estimator('fgo', make_config(window_size=12), data)
        self.assertGreater(np.max(np.abs(result['X']-result['X_online'])), .01)
        residual = measurement_residuals(result, data)
        self.assertTrue(np.isnan(residual[:, 0]).all())
        for k in range(1, 12):
            expected = [data['toa_measurements'][j,k]-range_measurement(result['X_online'][:,k],e)
                        for j,e in enumerate(data['emitter_positions'].T)]
            np.testing.assert_allclose(residual[:,k], expected)
        stats = metrics(result, data['true_positions'], elapsed)
        errors = result['X_online'][:2,1:]-data['true_positions'][:,1:]
        self.assertAlmostEqual(stats['rmse'], np.sqrt(np.mean(np.sum(errors**2,axis=0))))
        self.assertEqual(stats['processed_epochs'], 11)

    def test_residual_views_and_json_use_missing_not_zero(self):
        data = make_data(num_steps=8)
        result, elapsed = run_estimator('kfv', make_config(), data)
        item = dict(name='test', kind='kfv', result=result, metrics=metrics(result,data['true_positions'],elapsed))
        experiment = dict(data=data, estimators=[item,item])
        payload = comparison_payload(experiment)
        encoded = json.dumps(payload, allow_nan=False)
        for mode in ('postfit', 'prefit', 'truth'):
            self.assertIsNone(payload['residuals'][mode][0][0][0])
        self.assertNotEqual(payload['residuals']['prefit'], payload['residuals']['postfit'])
        self.assertNotIn('NaN',encoded)
        html = animate_pair(experiment).data
        self.assertNotIn('__PAYLOAD__',html)
        self.assertNotIn('__ID__',html)
        self.assertIn('Iteration view',html)

    def test_notebook_cells_compile_and_setup_is_non_destructive(self):
        notebook = json.loads(Path('KFV-FGO-Colab.ipynb').read_text())
        for cell in notebook['cells']:
            if cell['cell_type']=='code':
                compile(''.join(cell['source']), '<colab>', 'exec')
        setup=''.join(notebook['cells'][0]['source'])
        self.assertNotIn('rmtree',setup)
        self.assertIn('BRANCH = "main"',setup)

if __name__=='__main__':
    unittest.main()
