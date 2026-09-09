from __future__ import annotations

from dataclasses import dataclass, field
import numpy as np
from ..config.config import range_measurement, range_jacobian, motion, motion_jacobian


def robust_sqrt_weight(residual, loss_type, delta):
    r = abs(float(residual)); kind = loss_type.lower()
    if kind == "huber": value = 1.0 if r <= delta else np.sqrt(delta / r)
    elif kind == "cauchy": value = np.sqrt(1.0 / (1.0 + (r / delta) ** 2))
    elif kind == "tukey": value = np.sqrt((1.0 - (r / delta) ** 2) ** 2) if r <= delta else 0.0
    elif kind == "none": value = 1.0
    else: raise ValueError(f"Unknown loss type: {loss_type}")
    return value


@dataclass
class State:
    gid: int
    lid: int
    value: np.ndarray
    status: str = "Forward"


class Factor:
    def __init__(self, states, z=None, omega=None):
        self.states = list(states); self.z = z; self.omega = omega; self.A = None; self.b = None; self.status = ""
    def evaluate(self): raise NotImplementedError


def whiten(A, b, omega):
    L = np.linalg.cholesky(np.asarray(omega, dtype=float))
    return L.T @ A, L.T @ np.asarray(b).reshape(-1)


class PositionFactor(Factor):
    def evaluate(self):
        self.A, self.b = whiten(-np.eye(4), self.z - self.states[0].value, self.omega); return self


class RangeFactor(Factor):
    def evaluate(self):
        measurement, emitter = self.z["range"], self.z["emitter"]
        residual = measurement - range_measurement(self.states[0].value, emitter)
        jacobian = -range_jacobian(self.states[0].value, emitter)
        weight = robust_sqrt_weight(abs(residual) * np.sqrt(float(self.omega[0, 0])), self.z["loss_type"], self.z["loss_delta"])
        self.A, self.b = whiten(weight * jacobian[None, :], weight * residual, self.omega); return self


class PropagateFactor(Factor):
    def __init__(self, states, config): super().__init__(states); self.config = config
    def evaluate(self):
        cfg = self.config
        F = motion_jacobian(self.states[0].value, cfg.dt, cfg.omega)
        residual = self.states[1].value - motion(self.states[0].value, cfg.dt, cfg.omega)
        self.A, self.b = whiten(np.hstack((-F, np.eye(4))), residual, np.linalg.inv(cfg.q)); return self


class PdrFactor(Factor):
    def evaluate(self):
        residual = self.z - (self.states[1].value - self.states[0].value)
        self.A, self.b = whiten(np.hstack((np.eye(4), -np.eye(4))), residual, self.omega); return self


class MarginFactor(Factor):
    def __init__(self, states, A, b, omega=None):
        super().__init__(states, None, omega); self.A, self.b = A, np.asarray(b).reshape(-1)
    def evaluate(self): return self


class AutoDiffFactor(Factor):
    def __init__(self, states, z, omega, error_func, delta=1e-4):
        super().__init__(states, z, omega); self.error_func, self.delta = error_func, delta
    def evaluate(self):
        x = self.states[0].value; base = np.asarray(self.error_func(x, self.z), dtype=float).reshape(-1)
        jacobian = np.zeros((len(base), len(x)))
        for index in range(len(x)):
            plus, minus = x.copy(), x.copy(); plus[index] += self.delta; minus[index] -= self.delta
            jacobian[:, index] = (self.error_func(plus, self.z) - self.error_func(minus, self.z)) / (2 * self.delta)
        self.A, self.b = whiten(jacobian, base, self.omega); return self


class AutoDiffRobustFactor(AutoDiffFactor):
    def __init__(self, states, z, error_func, kernel_func):
        super().__init__(states, z, np.eye(1), error_func)
        self.kernel_func = kernel_func

    def evaluate(self):
        x = self.states[0].value
        residual = np.asarray(self.kernel_func(self.error_func(x, self.z))).reshape(-1)
        jacobian = np.zeros((len(residual), len(x)))
        for index in range(len(x)):
            plus, minus = x.copy(), x.copy(); plus[index] += self.delta; minus[index] -= self.delta
            jacobian[:, index] = (self.kernel_func(self.error_func(plus, self.z)) -
                                  self.kernel_func(self.error_func(minus, self.z))) / (2 * self.delta)
        self.A, self.b = jacobian, residual
        return self
