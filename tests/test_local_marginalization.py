"""Independent dense-Schur and posterior checks for local QR marginalization.

Run from the repository root: python -m unittest discover -s tests -v
"""
import unittest
from unittest.mock import patch

import numpy as np

from config.config import FgoConfig
from core.fgo.factor.factor import Factor, State
from core.fgo.factor.margin_factor import MarginFactor
from core.fgo.factor.position_factor import PositionFactor
from core.fgo.factor_graph import FactorGraph


class LinearFactor(Factor):
    def __init__(self, states, a, z):
        super().__init__(states, np.asarray(z, dtype=float))
        self.matrix = np.asarray(a, dtype=float)
        self.evaluations = 0

    def evaluate(self):
        self.evaluations += 1
        self.A = self.matrix
        self.b = self.z - self.A @ np.concatenate([s.value for s in self.states])
        return self


def chain(count=5):
    rng = np.random.default_rng(42)
    graph = FactorGraph(FgoConfig(max_iteration=10))
    for gid in range(1, count + 1):
        state = State(gid, gid, rng.normal(size=4))
        graph.add_state(state)
        graph.add_factor(LinearFactor([state], rng.normal(size=(6, 4)), rng.normal(size=6)))
        if gid == 1:
            graph.add_factor(LinearFactor([state], 2 * np.eye(4), rng.normal(size=4)))
        else:
            graph.add_factor(LinearFactor(graph.states[-2:],
                                          np.hstack((-3 * np.eye(4), 3 * np.eye(4))),
                                          rng.normal(size=4)))
    return graph


def dense_system(graph):
    """Assemble an independent reference using global state selection matrices."""
    states = graph.active_states
    eye = np.eye(4 * len(states))
    selections = {s.gid: eye[4 * i:4 * i + 4] for i, s in enumerate(states)}
    blocks, rhs = [], []
    for factor in graph.factors:
        if factor.status == "Margin":
            continue
        factor.evaluate()
        select = np.vstack([selections[s.gid] for s in factor.states])
        blocks.append(factor.A @ select)
        rhs.append(factor.b)
    return np.vstack(blocks), np.concatenate(rhs)


def dense_schur(j, r, removed_indices, rank_deficient=False):
    h, g = j.T @ j, j.T @ r
    retained = np.array([i for i in range(h.shape[0]) if i not in removed_indices])
    removed = np.array(removed_indices)
    hmm = h[np.ix_(removed, removed)]
    hmk = h[np.ix_(removed, retained)]
    if rank_deficient:
        solved_h, solved_g = np.linalg.pinv(hmm) @ hmk, np.linalg.pinv(hmm) @ g[removed]
    else:
        solved_h, solved_g = np.linalg.solve(hmm, hmk), np.linalg.solve(hmm, g[removed])
    return (h[np.ix_(retained, retained)] - hmk.T @ solved_h,
            g[retained] - hmk.T @ solved_g)


class LocalMarginalizationTests(unittest.TestCase):
    def assert_system(self, graph, expected_h, expected_g):
        graph.normal_equation()
        np.testing.assert_allclose(graph.J.T @ graph.J, expected_h, rtol=1e-11, atol=1e-11)
        np.testing.assert_allclose(graph.J.T @ graph.r, expected_g, rtol=1e-11, atol=1e-11)

    def test_matches_full_schur_without_absorbing_retained_factors(self):
        graph = chain()
        j, r = dense_system(graph)
        expected_h, expected_g = dense_schur(j, r, list(range(4)))
        unaffected = [f for f in graph.factors if all(s.gid != 1 for s in f.states)]
        for factor in graph.factors:
            factor.evaluations = 0
        graph.marginalize(1)
        self.assertTrue(all(f.evaluations == 0 for f in unaffected))
        self.assertTrue(all(sum(item is f for item in graph.factors) == 1 for f in unaffected))
        self.assert_system(graph, expected_h, expected_g)
        priors = [f for f in graph.factors if isinstance(f, MarginFactor)]
        self.assertEqual(len(priors), 1)
        self.assertEqual([s.gid for s in priors[0].states], [2])
        self.assertEqual(priors[0].A0.shape, (4, 4))

    def test_prior_remains_anchored_after_retained_states_move(self):
        graph = chain()
        j, r = dense_system(graph)
        h, g = dense_schur(j, r, list(range(4)))
        graph.marginalize(1)
        shift = np.linspace(-.7, .5, h.shape[0])
        for i, state in enumerate(graph.active_states):
            state.value += shift[4 * i:4 * i + 4]
        self.assert_system(graph, h, g - h @ shift)

    def test_repeated_removal_preserves_batch_posterior(self):
        graph = chain(7)
        j, r = dense_system(graph)
        full_delta = np.linalg.lstsq(j, r, rcond=None)[0]
        for gid in (1, 2, 3, 4):
            graph.marginalize(gid)
            h, g = dense_schur(j, r, list(range(4 * gid)))
            self.assert_system(graph, h, g)
            remaining_delta = np.linalg.lstsq(graph.J, graph.r, rcond=None)[0]
            np.testing.assert_allclose(remaining_delta, full_delta[4 * gid:], atol=1e-11)
            self.assertEqual(sum(isinstance(f, MarginFactor) for f in graph.factors), 1)

    def test_multiple_and_nonconsecutive_removed_states(self):
        for gids in ([1, 2], [2, 4]):
            with self.subTest(gids=gids):
                graph = chain(6)
                j, r = dense_system(graph)
                removed = [4 * (gid - 1) + offset for gid in gids for offset in range(4)]
                h, g = dense_schur(j, r, removed)
                graph.marginalize(gids)
                self.assert_system(graph, h, g)

    def test_rank_deficient_elimination_matches_generalized_schur(self):
        graph = FactorGraph(FgoConfig())
        states = [State(i, i, np.arange(4, dtype=float) * i) for i in (1, 2)]
        for state in states:
            graph.add_state(state)
        partial = np.eye(4)[:2]
        graph.add_factor(LinearFactor([states[0]], 2 * partial, [1, 3]))
        graph.add_factor(LinearFactor(states, np.hstack((-partial, partial)), [2, -1]))
        graph.add_factor(LinearFactor([states[1]], np.eye(4), [1, 2, 3, 4]))
        j, r = dense_system(graph)
        h, g = dense_schur(j, r, list(range(4)), rank_deficient=True)
        graph.marginalize(1)
        self.assert_system(graph, h, g)
        self.assertEqual(graph.last_marginalization["removed_rank"], 2)

    def test_no_global_assembly_or_svd_and_constant_boundary_size(self):
        for window in (3, 40):
            graph = chain(window)
            with patch.object(graph, "normal_equation", side_effect=AssertionError("global assembly")), \
                    patch.object(np.linalg, "svd", side_effect=AssertionError("SVD")):
                graph.marginalize(1)
                graph.marginalize(2)
            self.assertEqual(graph.last_marginalization["boundary_dim"], 4)
            self.assertEqual(graph.last_marginalization["local_rows"], 14)
            self.assertEqual(graph.last_marginalization["prior_rows"], 4)
            self.assertEqual(len(graph.states), window)
            self.assertEqual(graph.win_size, window - 2)
            self.assertTrue(all(f.status != "Margin" for f in graph.factors))
            self.assertTrue(all(s.status != "Margin" for f in graph.factors for s in f.states))

    def test_discard_drops_incident_factors_without_a_prior(self):
        graph = chain(5)
        expected = [f for f in graph.factors if all(s.gid != 1 for s in f.states)]
        with patch.object(LinearFactor, "evaluate", side_effect=AssertionError("factor evaluated")):
            graph.discard(1)
        self.assertEqual(graph.factors, expected)
        self.assertEqual([s.lid for s in graph.active_states], [1, 2, 3, 4])
        self.assertEqual(len(graph.states), 5)
        self.assertEqual(graph.states[0].status, "Margin")

    def test_empty_removal_and_removal_without_boundary(self):
        graph = chain(2)
        factors = list(graph.factors)
        graph.marginalize(99)
        self.assertEqual(graph.factors, factors)
        graph.marginalize([1, 2])
        self.assertEqual(graph.factors, [])
        self.assertEqual(graph.active_states, [])
        self.assertEqual(len(graph.states), 2)

    def test_position_prior_moves_state_towards_measurement(self):
        graph = FactorGraph(FgoConfig(max_iteration=10))
        state = State(1, 1, np.array([10., -5., 2., 4.]))
        target = np.array([1., 2., 3., 4.])
        graph.add_state(state)
        graph.add_factor(PositionFactor([state], target, np.diag([1., 2., 3., 4.])))
        graph.estimate()
        np.testing.assert_allclose(state.value, target, atol=1e-12)


class EstimatorIntegrationTests(unittest.TestCase):
    def test_window_retirement_and_batch_policy_identity(self):
        from schur_window_benchmark import run_case
        from data.circle_eval import generate_data
        data = generate_data(num_steps=8, gmm_weights=(1., 0.), seed=7)
        trajectories = []
        for policy in ("schur", "discard"):
            record, trajectory = run_case(data, 8, policy, profile_stages=True)
            self.assertEqual(record["removal_calls"], 0)
            self.assertEqual(record["final_active_states"], 8)
            self.assertLessEqual(record["total_gn_iterations"], 7 * 10)
            self.assertEqual(record["marginalize_ms"], 0)
            for stage in ("add_state_factor_ms", "estimate_ms", "marginalize_ms"):
                self.assertAlmostEqual(sum(row[stage] for row in record["_stage_rows"]), record[stage])
            trajectories.append(trajectory)
            small, _ = run_case(data, 2, policy)
            self.assertEqual(small["removal_calls"], 6)
            self.assertEqual(small["final_active_states"], 2)
            self.assertEqual(small["max_solve_states"], 3)
        np.testing.assert_array_equal(*trajectories)


if __name__ == "__main__":
    unittest.main()
