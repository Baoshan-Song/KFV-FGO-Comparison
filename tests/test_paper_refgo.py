"""Independent algebraic checks of the v2 prediction, update and QR prior."""
import importlib
import unittest
from unittest.mock import patch
import numpy as np
from scipy.linalg import qr
from config.config import (Config, FilterConfig, FgoConfig, convert_kfv_to_fgo,
                           motion, motion_jacobian)
from core.estimator import KfvEstimator, FgoEstimator
from core.fgo.factor import Factor, State, MarginFactor, RangeFactor, PositionFactor, PropagateFactor
from core.fgo.factor_graph import FactorGraph
from core.filter.miekf import miekf
from data.circle_eval import generate_data


class LinearObservation(Factor):
    def evaluate(self):
        self.A = np.array([[1., 0., 0., 0.]])
        self.b = np.array([1.-self.states[0].value[0]])
        return self


class PaperConsistencyTests(unittest.TestCase):
    def test_fixed_prediction_prevents_repeated_measurement_counting(self):
        # Exact posterior: first component mean .5, variance .5, for any
        # iteration budget. A frozen prior RHS incorrectly tends to mean 1.
        for iterations in (1, 2, 10):
            with self.subTest(iterations=iterations):
                graph = FactorGraph(FgoConfig(max_iteration=iterations, threshold_iteration=0.))
                state = State(1, 1, np.zeros(4))
                graph.add_state(state).add_factor(MarginFactor([state], np.eye(4), np.zeros(4), np.zeros(4)))
                graph.add_factor(LinearObservation([state])).estimate().mar_measurements(1)
                self.assertAlmostEqual(state.value[0], .5, places=13)
                self.assertAlmostEqual(graph.last_measurement_marginalization['covariance'][0, 0], .5, places=13)
                x, p, *_ = miekf(np.zeros(4), np.eye(4), 1., 0.,
                    lambda x, dt, w: x.copy(), lambda *args: np.eye(4), np.zeros((4, 4)),
                    np.ones(1), np.zeros((2, 1)), lambda x, e: x[0],
                    lambda x, e: np.array([1., 0., 0., 0.]), 1., iterations, 0.)
                self.assertAlmostEqual(x[0], .5, places=13)
                self.assertAlmostEqual(p[0, 0], .5, places=13)

    def test_stage_one_is_covariance_prediction_before_observations(self):
        rng = np.random.default_rng(13)
        for _ in range(5):
            mat = rng.normal(size=(4, 4)); covariance = mat@mat.T + np.eye(4)
            cfg = FgoConfig()
            old = State(1, 1, rng.normal(size=4)); new = State(2, 2, motion(old.value, cfg.dt, cfg.omega))
            graph = FactorGraph(cfg)
            graph.add_state(old).add_factor(PositionFactor([old], old.value.copy(), np.linalg.inv(covariance)))
            graph.add_state(new).add_factor(PropagateFactor([old, new], cfg)).marginalize(1)
            prior = graph.factors[0].evaluate(); F = motion_jacobian(old.value, cfg.dt, cfg.omega)
            expected = np.linalg.inv(F@covariance@F.T + cfg.q)
            np.testing.assert_allclose(prior.A.T@prior.A, expected, atol=1e-12)
            np.testing.assert_allclose(prior.b, 0., atol=1e-12)
            self.assertEqual(len(graph.active_states), 1)

    def test_second_stage_uses_qr_preserves_affine_rhs_and_full_cost(self):
        graph = FactorGraph(FgoConfig(max_iteration=1))
        state = State(1, 1, np.zeros(4))
        graph.add_state(state).add_factor(MarginFactor([state], np.eye(4), np.zeros(4), np.zeros(4)))
        graph.add_factor(LinearObservation([state])).estimate()
        system = graph.last_linear_system
        # Move the state before compression: QR must preserve the actual
        # linear system, not manufacture a zero residual at this arbitrary mean.
        state.value += .2
        with patch('core.fgo.factor_graph.qr', wraps=qr) as decomposition:
            graph.mar_measurements(1)
            decomposition.assert_called_once()
        prior = graph.factors[0]; info = graph.last_measurement_marginalization
        np.testing.assert_allclose(prior.A0.T@prior.A0, system['A'].T@system['A'], atol=1e-13)
        np.testing.assert_allclose(prior.A0.T@prior.b0, system['A'].T@system['b'], atol=1e-13)
        self.assertGreater(np.linalg.norm(prior.b0), .1)
        self.assertGreater(np.linalg.norm(prior.evaluate().b), .1)
        for delta in np.random.default_rng(0).normal(size=(10, 4)):
            full = .5*np.linalg.norm(system['A']@delta-system['b'])**2
            compressed = .5*np.linalg.norm(prior.A0@delta-prior.b0)**2 + info['constant_cost']
            self.assertAlmostEqual(full, compressed, places=12)
        self.assertIs(graph.active_states[0], state)
        self.assertEqual(len(graph.factors), 1)

    def test_robust_weights_never_modify_nominal_noise_or_accumulate(self):
        state = State(1, 1, np.array([1., 0., 0., 0.]))
        for kernel in ('huber', 'cauchy', 'tukey', 'none'):
            omega = np.array([[100.]])
            factor = RangeFactor([state], dict(range=3., emitter=np.zeros(2), loss_type=kernel, loss_delta=2.), omega)
            factor.evaluate(); a, b = factor.A.copy(), factor.b.copy()
            for _ in range(4):
                factor.evaluate()
                np.testing.assert_array_equal(factor.A, a)
                np.testing.assert_array_equal(factor.b, b)
            np.testing.assert_array_equal(omega, [[100.]])

    def test_four_paper_scenes_match_independent_gain_filter(self):
        for radius in (105., 1000.):
            for outliers in (0., .2):
                data = generate_data(emitter_radius=radius, gmm_weights=(1-outliers, outliers), seed=7)
                for mode in ('EKF', 'iEKF', 'rEKF', 'riEKF'):
                    with self.subTest(radius=radius, outliers=outliers, mode=mode):
                        config = convert_kfv_to_fgo(Config(kfv=FilterConfig(mode=mode, max_iteration=10)))
                        result = KfvEstimator(config, data).run()
                        with patch.object(KfvEstimator, 'run', side_effect=AssertionError('Graph called filter')):
                            graph = FgoEstimator(config, data).run()
                        np.testing.assert_allclose(graph['X'], result['X'], rtol=0, atol=1e-8)
                        np.testing.assert_allclose(graph['covariances'], result['covariances'], rtol=1e-8, atol=1e-10)
                        for f, k in zip(graph['epoch_diagnostics'][1:], result['epoch_diagnostics'][1:]):
                            self.assertEqual(f['iterations'], k['iterations'])
                            np.testing.assert_allclose([t['objective'] for t in f['iteration_trace']],
                                [t['objective'] for t in k['iteration_trace']], rtol=1e-8, atol=1e-7)

    def test_alternative_kernels_and_iteration_budgets(self):
        data = generate_data(num_steps=25, seed=11, gmm_weights=(.8,.2))
        for kernel in ('huber', 'cauchy', 'tukey'):
            for iterations in (1, 30):
                config = convert_kfv_to_fgo(Config(kfv=FilterConfig(mode='riEKF',
                    max_iteration=iterations, robust_kernel=kernel, err_x0=np.array([1., -1., 0., 0.]))))
                a, b = KfvEstimator(config, data).run(), FgoEstimator(config, data).run()
                np.testing.assert_allclose(a['X'], b['X'], atol=1e-8, rtol=0)
                # Tukey can reject every observation, leaving a large prediction
                # covariance. Use a matrix-scaled tolerance, including near-zero entries.
                for pa, pb in zip(a['covariances'], b['covariances']):
                    self.assertLess(np.linalg.norm(pa-pb), 1e-10*max(1., np.linalg.norm(pa)))

    def test_one_iteration_does_not_relinearize_at_posterior(self):
        config = convert_kfv_to_fgo(Config(kfv=FilterConfig(mode='EKF')))
        result = FgoEstimator(config, generate_data(num_steps=2)).run()
        self.assertEqual(result['epoch_diagnostics'][1]['iterations'], 1)
        self.assertFalse(result['epoch_diagnostics'][1]['converged'])
        # Large first correction makes an unauthorized second linearization visible.
        expected = KfvEstimator(config, generate_data(num_steps=2)).run()
        np.testing.assert_allclose(result['covariances'], expected['covariances'], atol=1e-11)

    def test_numerical_jacobian_switch_is_used(self):
        config = Config(fgo=FgoConfig(autodiff=True, max_iteration=1))
        with patch('core.fgo.factor.RangeFactor.range_jacobian', side_effect=AssertionError('Analytic Jacobian called')):
            FgoEstimator(config, generate_data(num_steps=3)).run()

    def test_refgo_rejects_larger_window_and_no_solve_qr(self):
        with self.assertRaisesRegex(ValueError, 'window_size=1'):
            FgoEstimator(Config(fgo=FgoConfig(window_size=2)), generate_data(num_steps=3)).run()
        graph = FactorGraph(FgoConfig()); graph.add_state(State(1, 1, np.zeros(4)))
        with self.assertRaisesRegex(ValueError, 'Estimate'):
            graph.mar_measurements(1)

    def test_public_api_cannot_dispatch_fgo_to_filter(self):
        import kfv_fgo
        from kfv_fgo.estimators import FgoEstimator as public_fgo
        self.assertIs(public_fgo, FgoEstimator)
        self.assertIs(kfv_fgo.FgoEstimator, FgoEstimator)
        with patch.object(KfvEstimator, 'run', side_effect=AssertionError('Graph called filter')):
            public_fgo(Config(), generate_data(num_steps=3)).run()

if __name__ == '__main__':
    unittest.main()
