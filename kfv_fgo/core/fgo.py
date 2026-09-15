"""Factor signs and frozen marginalization priors follow the MATLAB reference.

The graph solves A * delta = b and ADDS delta. A is not uniformly the
residual derivative in this research implementation; do not change one sign alone.
"""

from dataclasses import dataclass
from time import perf_counter

import numpy as np

from ..model.simulation import (
    motion,
    motion_jacobian,
    range_jacobian,
    range_measurement,
)
from .filter import _weight


@dataclass
class State:
    gid: int
    lid: int
    value: np.ndarray
    status: str = "Forward"


class Factor:
    def __init__(self, states, z=None, omega=None):
        self.states = list(states)
        self.z, self.omega = z, omega
        self.A = self.b = None
        self.status = ""

    def evaluate(self):
        raise NotImplementedError


def whiten(A, b, omega):
    root = np.linalg.cholesky(np.asarray(omega)).T
    return root @ A, root @ np.asarray(b).reshape(-1)


def robust_sqrt_weight(residual, loss_type, delta):
    return np.sqrt(_weight(residual, loss_type, delta))


class PositionFactor(Factor):
    def evaluate(self):
        self.A, self.b = whiten(
            -np.eye(len(self.states[0].value)),
            self.z - self.states[0].value,
            self.omega,
        )
        return self


class RangeFactor(Factor):
    def evaluate(self):
        state = self.states[0].value
        residual = self.z["range"] - range_measurement(state, self.z["emitter"])
        information = float(np.asarray(self.omega).reshape(-1)[0])
        weight = _weight(
            abs(residual) * np.sqrt(information),
            self.z["loss_type"],
            self.z["loss_delta"],
        )
        root = np.sqrt(information * weight)
        self.A = root * range_jacobian(state, self.z["emitter"])[None, :]
        self.b = np.array([root * residual])
        # MATLAB RangeFactor mutates Omega on every evaluation, even marginalization.
        self.omega = np.array([[information * weight]])
        return self


class PropagateFactor(Factor):
    def __init__(self, states, config, propagate=None):
        super().__init__(states)
        self.config, self.propagate = config, propagate

    def evaluate(self):
        old, new = self.states
        if self.propagate is None:
            cfg = self.config
            predicted = motion(old.value, cfg.dt, cfg.omega)
            F, Q = motion_jacobian(old.value, cfg.dt, cfg.omega), cfg.q
        else:
            predicted, F, Q = self.propagate(old.value)
        self.A, self.b = whiten(
            np.hstack((F, -np.eye(len(new.value)))),
            new.value - predicted,
            np.linalg.inv(Q),
        )
        return self


class GnssPseudorangeFactor(Factor):
    def __init__(self, states, observe, config):
        super().__init__(states)
        self.observe, self.config = observe, config

    def evaluate(self):
        residual, H, R = self.observe(self.states[0].value)[:3]
        weights = np.array(
            [
                _weight(v, self.config.robust_kernel, self.config.robust_delta)
                for v in np.abs(residual) / np.sqrt(np.diag(R))
            ]
        )
        root = np.sqrt(weights)
        lower = np.linalg.cholesky(R)
        self.A = root[:, None] * np.linalg.solve(lower, H)
        self.b = root * np.linalg.solve(lower, residual)
        self.omega = np.diag(weights) @ np.linalg.inv(R)
        return self


class PdrFactor(Factor):
    def evaluate(self):
        n = len(self.states[0].value)
        self.A, self.b = whiten(
            np.hstack((-np.eye(n), np.eye(n))),
            self.z - (self.states[1].value - self.states[0].value),
            self.omega,
        )
        return self


class MarginFactor(Factor):
    """Frozen A and b intentionally reproduce MATLAB, not a reanchored prior."""

    def __init__(self, states, A, b, omega=None):
        super().__init__(states, None, omega)
        self.A, self.b = A.copy(), np.asarray(b).reshape(-1).copy()

    def evaluate(self):
        return self


class AutoDiffFactor(Factor):
    """Legacy name: central finite differences, not an autodiff engine."""

    def __init__(self, states, z, omega, error_func, delta=1e-4):
        super().__init__(states, z, omega)
        self.error_func, self.delta = error_func, delta

    def evaluate(self):
        x = self.states[0].value
        residual = np.asarray(self.error_func(x, self.z)).reshape(-1)
        jacobian = np.empty((len(residual), len(x)))
        for i in range(len(x)):
            plus, minus = x.copy(), x.copy()
            plus[i] += self.delta
            minus[i] -= self.delta
            jacobian[:, i] = (
                self.error_func(plus, self.z) - self.error_func(minus, self.z)
            ) / (2 * self.delta)
        self.A, self.b = whiten(-jacobian, residual, self.omega)
        return self


class AutoDiffRobustFactor(AutoDiffFactor):
    def __init__(self, states, z, error_func, kernel_func):
        super().__init__(
            states, z, np.eye(1), lambda x, z: kernel_func(error_func(x, z))
        )

    def evaluate(self):
        count = np.asarray(self.error_func(self.states[0].value, self.z)).size
        self.omega = np.eye(count)
        return super().evaluate()


class FactorGraph:
    def __init__(self, config, trace=False):
        self.config, self.trace = config, trace
        self.states, self.factors = [], []
        self.J = self.r = self.latest_information_matrix = None
        self.ls_time = self.margin_time = self.margin_meas_time = 0.0
        self.residual_norm_all, self.iterations = [], []

    @property
    def active_states(self):
        return [s for s in self.states if s.status != "Margin"]

    @property
    def win_size(self):
        return len(self.active_states)

    def add_state(self, state):
        state.lid = self.win_size + 1
        self.states.append(state)
        return self

    def add_factor(self, factor):
        self.factors.append(factor)
        return self

    def normal_equation(self):
        factors = [f.evaluate() for f in self.factors if f.status != "Margin"]
        states = self.active_states
        n = len(states[0].value)
        offsets = {s.gid: i * n for i, s in enumerate(states)}
        self.J = np.zeros((sum(len(f.b) for f in factors), n * len(states)))
        self.r = np.zeros(self.J.shape[0])
        row = 0
        for f in factors:
            rows = slice(row, row + len(f.b))
            self.r[rows] = f.b
            for j, s in enumerate(f.states):
                col = offsets[s.gid]
                self.J[rows, col : col + n] = f.A[:, j * n : (j + 1) * n]
            row += len(f.b)
        return self

    def estimate(self):
        self.iterations, self.residual_norm_all = [], []
        started = perf_counter()
        for iteration in range(self.config.max_iteration):
            self.normal_equation()
            information = self.J.T @ self.J
            delta = np.linalg.solve(information, self.J.T @ self.r)
            states = self.active_states
            self.residual_norm_all.append(
                np.r_[states[0].value.copy(), np.linalg.norm(self.r)]
            )
            if self.trace:
                self.iterations.append(
                    {
                        "J": self.J.copy(),
                        "r": self.r.copy(),
                        "delta": delta.copy(),
                        "information": information.copy(),
                    }
                )
            n = len(states[0].value)
            for i, s in enumerate(states):
                s.value = s.value + delta[i * n : (i + 1) * n]
                if not np.all(np.isfinite(s.value)):
                    raise FloatingPointError("FGO produced non-finite state")
            self.latest_information_matrix = information
            if np.linalg.norm(delta) / len(delta) < self.config.threshold_iteration:
                break
        self.iteration_count = iteration + 1
        self.ls_time = perf_counter() - started
        return self

    def marginalize(self, gids):
        started = perf_counter()
        remove = {gids} if np.isscalar(gids) else set(gids)
        self.normal_equation()
        states = self.active_states
        n = len(states[0].value)
        removed = [i for i, s in enumerate(states) if s.gid in remove]
        remaining = [i for i, s in enumerate(states) if s.gid not in remove]
        if not removed or not remaining:
            raise ValueError("Marginalization needs removed and remaining states")
        c1 = np.concatenate([np.arange(i * n, (i + 1) * n) for i in removed])
        c2 = np.concatenate([np.arange(i * n, (i + 1) * n) for i in remaining])
        J1, J2 = self.J[:, c1], self.J[:, c2]
        H11, H12, H21, H22 = J1.T @ J1, J1.T @ J2, J2.T @ J1, J2.T @ J2
        np.linalg.cholesky(H11)
        h = H22 - H21 @ np.linalg.solve(H11, H12)
        b = J2.T @ self.r - H21 @ np.linalg.solve(H11, J1.T @ self.r)
        U, singular, _ = np.linalg.svd(h)
        root = np.sqrt(singular)[:, None] * U.T
        prior = np.linalg.lstsq(root.T, b, rcond=None)[0]
        for s in states:
            if s.gid in remove:
                s.status, s.lid = "Margin", 0
        # Intentionally retain unaffected factors, as MATLAB does, even though
        # their information also participates in the Schur prior. See docs.
        self.factors = [
            f for f in self.factors if not any(s.gid in remove for s in f.states)
        ]
        for i, s in enumerate(self.active_states, 1):
            s.lid = i
        self.add_factor(MarginFactor([states[i] for i in remaining], root, prior))
        self.margin_time = perf_counter() - started
        return self

    def mar_measurements(self, gid):
        started = perf_counter()
        state = next(s for s in self.states if s.gid == gid)
        self.factors = [
            f for f in self.factors if not any(s.gid == gid for s in f.states)
        ]
        # Reuse the LAST solve's information: do not relinearize here.
        self.add_factor(
            PositionFactor(
                [state], state.value.copy(), self.latest_information_matrix.copy()
            )
        )
        self.margin_meas_time = perf_counter() - started
        return self
