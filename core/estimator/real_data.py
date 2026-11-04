"""GNSS/INS adapters for independent KFV and square-root FGO solvers.

The observation and first-sample Euler propagation models are from ZzPolyU's
PR #5. Re-FGO uses both canonical QR marginalizations, never a filter call.
"""
from time import perf_counter
import numpy as np
from config.config import validate_config, _config_to_dict
from core.filter.miekf import observation_update
from core.fgo.factor_graph import FactorGraph
from core.fgo.factor import State, PositionFactor
from core.fgo.factor.factor import Factor
from core.objective import robust_weights, loss_cost


class ImuFactor(Factor):
    def __init__(self, states, propagate):
        super().__init__(states)
        self.propagate = propagate

    def evaluate(self):
        prediction, jacobian, covariance = self.propagate(self.states[0].value)
        root = np.linalg.cholesky(covariance)
        self.A = np.linalg.solve(root, np.column_stack((-jacobian, np.eye(len(prediction)))))
        self.b = np.linalg.solve(root, prediction-self.states[1].value)
        return self


class GnssFactor(Factor):
    def __init__(self, state, observe, cfg):
        super().__init__([state])
        self.observe, self.cfg = observe, cfg

    def evaluate(self):
        raw, H, covariance, self.ids = self.observe(self.states[0].value)
        sigma = np.sqrt(np.diag(covariance))
        weights = robust_weights(raw/sigma, self.cfg.robust_kernel, self.cfg.robust_delta)
        scale = np.sqrt(weights)/sigma
        self.A, self.b = scale[:, None]*H, scale*raw
        return self

    def observation_diagnostic(self):
        raw, _, covariance, _ = self.observe(self.states[0].value)
        return loss_cost(raw/np.sqrt(np.diag(covariance)),
                         self.cfg.robust_kernel, self.cfg.robust_delta), raw


def allocate(config, data, cfg):
    x0 = data.initial_state + cfg.err_x0
    validate_config(config, x0, cfg)
    n, d = data.num_steps, len(x0)
    online = np.full((d, n), np.nan)
    predicted = np.full_like(online, np.nan)
    covariances = np.full((n, d, d), np.nan)
    online[:, 0], covariances[0] = x0, cfg.p0
    return x0, online, predicted, covariances, [None]*n


def finish(config, data, solver, start, online, predicted, covariance, diagnostics, final=None):
    return dict(X=online.copy() if final is None else final, X_online=online,
        X_pred=predicted, covariances=covariance, epoch_diagnostics=diagnostics,
        debug_info=diagnostics, solver=solver, timestamps=data.timestamps.copy(),
        runtime_seconds=perf_counter()-start, config=_config_to_dict(config),
        backend=data.backend.metadata,
        trajectory_semantics='X_online is causal; X is the final retained graph history')


def run_kfv(config, data):
    start = perf_counter()
    cfg = config.kfv
    state, online, predicted, covariances, diagnostics = allocate(config, data, cfg)
    covariance = cfg.p0.copy()
    if cfg.mode not in ('EKF', 'iEKF', 'rEKF', 'riEKF'):
        raise ValueError(f'Unknown filter {cfg.mode}')
    iterations = cfg.max_iteration if cfg.mode in ('iEKF', 'riEKF') else 1
    loss = cfg.robust_kernel if cfg.mode in ('rEKF', 'riEKF') else 'none'
    for k in range(1, data.num_steps):
        xp, F, Q = data.propagate(state, k, cfg)
        pp = F@covariance@F.T+Q
        predicted[:, k] = xp
        state, covariance, diagnostics[k] = observation_update(
            xp, pp, lambda x: data.linearize(x, k), iterations,
            cfg.threshold_iteration, loss, cfg.robust_delta)
        online[:, k], covariances[k] = state, covariance
    return finish(config, data, 'KfvEstimator', start, online, predicted, covariances, diagnostics)


def run_fgo(config, data):
    start = perf_counter()
    cfg = config.fgo
    x0, online, predicted, covariances, diagnostics = allocate(config, data, cfg)
    recursive = cfg.imitate_kfv
    if recursive and cfg.window_size != 1:
        raise ValueError('Re-FGO requires window_size=1')
    graph = FactorGraph(cfg)
    state = State(1, 1, x0.copy())
    graph.add_state(state)
    graph.add_factor(PositionFactor([state], x0.copy(), np.linalg.inv(cfg.p0)))
    for k in range(1, data.num_steps):
        old = graph.states[-1]
        xp, _, _ = data.propagate(old.value, k, cfg)
        predicted[:, k] = xp
        state = State(k+1, graph.win_size+1, xp.copy())
        graph.add_state(state)
        graph.add_factor(ImuFactor([old, state], lambda x, i=k: data.propagate(x, i, cfg)))
        if recursive:
            graph.marginalize(k)
        graph.add_factor(GnssFactor(state, lambda x, i=k: data.linearize(x, i), cfg))
        graph.estimate()
        online[:, k] = state.value
        raw, _, _, ids = data.linearize(state.value, k)
        prefit, _, _, pre_ids = data.linearize(xp, k)
        diagnostics[k] = dict(iteration_trace=graph.iteration_trace,
            iterations=len(graph.iteration_trace)-1, converged=graph.converged,
            postfit_residual=raw, residual=raw, postfit_ids=ids,
            prefit_residual=prefit, prefit_ids=pre_ids,
            active_states=graph.win_size)
        if recursive:
            graph.mar_measurements(k+1)
            covariances[k] = graph.last_measurement_marginalization['covariance']
        elif graph.win_size > cfg.window_size:
            graph.marginalize(k+1-cfg.window_size)
    final = np.column_stack([state.value for state in graph.states])
    return finish(config, data, 'FgoEstimator', start, online, predicted,
                  covariances, diagnostics, final)
