"""MATLAB recurrence compatibility, deliberately distinct from fixed-prior MAP."""
import importlib
import unittest
from unittest.mock import patch

import numpy as np

from config.config import Config, FgoConfig, FilterConfig
from core.estimator.FgoEstimator import FgoEstimator
from core.estimator.KfvEstimator import KfvEstimator
from core.fgo.factor.factor import Factor, State
from core.fgo.factor.margin_factor import MarginFactor
from core.fgo.factor.RangeFactor import RangeFactor
from core.fgo.factor_graph import FactorGraph
from data.circle_eval import generate_data


class LinearObservation(Factor):
    def evaluate(self):
        self.A = np.array([[1., 0., 0., 0.]])
        self.b = np.array([1. - self.states[0].value[0]])
        return self


class MatlabCompatibilityTests(unittest.TestCase):
    def test_legacy_recurrence_is_distinct_from_fixed_prior_map(self):
        # Prior N(0,I), observation x[0]=1 with variance 1: MAP is x[0]=0.5.
        # MATLAB's frozen prior residual instead produces 1 - 0.5**iterations.
        for compat, expected in ((False, .5), (True, 1. - .5**10)):
            with self.subTest(matlab_compat=compat):
                graph = FactorGraph(FgoConfig(max_iteration=10, threshold_iteration=0.))
                state = State(1, 1, np.zeros(4))
                graph.add_state(state)
                graph.add_factor(MarginFactor([state], np.eye(4), np.zeros(4),
                                              np.zeros(4), matlab_compat=compat))
                graph.add_factor(LinearObservation([state]))
                graph.estimate()
                self.assertAlmostEqual(state.value[0], expected, places=13)

    def test_legacy_weight_history_does_not_change_other_factors_or_base_noise(self):
        state = State(1, 1, np.array([1., 0., 0., 0.]))
        measurement = {"range": 3., "emitter": np.zeros(2),
                       "loss_type": "huber", "loss_delta": 2.}
        base_information = np.array([[100.]])
        legacy = RangeFactor([state], measurement, base_information, matlab_compat=True)
        standard = RangeFactor([state], measurement, base_information)
        standard.evaluate()
        initial_a, initial_b = standard.A.copy(), standard.b.copy()
        legacy.evaluate()
        np.testing.assert_allclose(legacy.A, initial_a)
        first_information = legacy.omega.copy()
        legacy.evaluate()
        self.assertLess(legacy.omega[0, 0], first_information[0, 0])
        np.testing.assert_array_equal(base_information, [[100.]])
        standard.evaluate()
        np.testing.assert_array_equal(standard.A, initial_a)
        np.testing.assert_array_equal(standard.b, initial_b)
        new_epoch = RangeFactor([state], measurement, base_information, matlab_compat=True)
        new_epoch.evaluate()
        np.testing.assert_array_equal(new_epoch.A, initial_a)

    def test_refgo_matches_independent_riekf_states_and_covariances(self):
        kfv_module = importlib.import_module("core.estimator.KfvEstimator")
        fgo_module = importlib.import_module("core.estimator.FgoEstimator")
        filter_update = kfv_module.rmiekf
        scenes = (("clean", 0., 0., 10.), ("entry", .2, 0., 10.),
                  ("notebook", .35, 30., 5.))
        for scene, outlier_weight, mean, sigma in scenes:
            data = generate_data(num_steps=100, seed=7,
                                 gmm_weights=(1. - outlier_weight, outlier_weight),
                                 gmm_means=(0., mean), gmm_sigmas=(.1, sigma))
            for iterations in (1, 10):
                with self.subTest(scene=scene, iterations=iterations):
                    config = Config(
                        kfv=FilterConfig(mode="riEKF", max_iteration=iterations,
                                         robust_kernel="huber", r=.01),
                        fgo=FgoConfig(imitate_kfv=True, window_size=1,
                                      max_iteration=iterations, robust_kernel="huber", r=.01))
                    kfv_covariances, fgo_covariances = [], []

                    def record_update(*args, **kwargs):
                        result = filter_update(*args, **kwargs)
                        kfv_covariances.append(result[1].copy())
                        return result

                    class RecordingGraph(FactorGraph):
                        def mar_measurements(self, gid):
                            fgo_covariances.append(np.linalg.inv(self.latest_information_matrix))
                            return super().mar_measurements(gid)

                    with patch.object(kfv_module, "rmiekf", side_effect=record_update):
                        kfv = KfvEstimator(config, data).run()["X"]
                    with patch.object(fgo_module, "FactorGraph", RecordingGraph), \
                            patch.object(kfv_module, "rmiekf",
                                         side_effect=AssertionError("FGO called the filter")):
                        fgo = FgoEstimator(config, data).run()["X"]
                    self.assertEqual(len(fgo_covariances), 99)
                    np.testing.assert_allclose(fgo, kfv, rtol=0., atol=1e-8)
                    np.testing.assert_allclose(fgo_covariances, kfv_covariances,
                                               rtol=1e-8, atol=1e-9)

    def test_standard_sw_fgo_does_not_enable_matlab_factors(self):
        fgo_module = importlib.import_module("core.estimator.FgoEstimator")
        observed_graphs = []

        class RecordingGraph(FactorGraph):
            def __init__(self, config):
                super().__init__(config)
                observed_graphs.append(self)

            def normal_equation(self):
                for factor in self.factors:
                    if isinstance(factor, (RangeFactor, MarginFactor)):
                        self.assert_standard(factor)
                return super().normal_equation()

            @staticmethod
            def assert_standard(factor):
                if factor.matlab_compat:
                    raise AssertionError("MATLAB factor leaked into standard SW-FGO")

        config = Config(fgo=FgoConfig(imitate_kfv=False, window_size=2,
                                      max_iteration=10, robust_kernel="huber"))
        data = generate_data(num_steps=12, seed=7)
        with patch.object(fgo_module, "FactorGraph", RecordingGraph):
            result = FgoEstimator(config, data).run()
        self.assertTrue(np.isfinite(result["X"]).all())
        self.assertEqual(observed_graphs[0].win_size, 2)
        for factor in observed_graphs[0].factors:
            if isinstance(factor, RangeFactor):
                np.testing.assert_array_equal(factor.omega, [[1. / config.fgo.r]])
            elif isinstance(factor, MarginFactor):
                self.assertFalse(factor.matlab_compat)

    def test_refgo_rejects_larger_windows(self):
        config = Config(fgo=FgoConfig(imitate_kfv=True, window_size=2))
        with self.assertRaisesRegex(ValueError, "window_size=1"):
            FgoEstimator(config, generate_data(num_steps=3)).run()

    def test_comparison_entry_uses_requested_configuration(self):
        import example_kfv_fgo_comparison as entry
        with patch.object(entry, "run_and_display", return_value=None) as display:
            entry.main()
        data, first, second = display.call_args.args
        self.assertEqual(data["num_steps"], 100)
        self.assertIs(first["config"], second["config"])
        config = first["config"]
        self.assertEqual(config.kfv.mode, "riEKF")
        self.assertTrue(config.fgo.imitate_kfv)
        for cfg in (config.kfv, config.fgo):
            self.assertEqual(cfg.max_iteration, 10)
            self.assertEqual(cfg.window_size, 1)
            self.assertEqual(cfg.robust_kernel, "huber")
            self.assertAlmostEqual(cfg.r, .01)


if __name__ == "__main__":
    unittest.main()
