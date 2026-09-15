"""Four-state simulation propagation and range observations."""

import numpy as np


def motion(x: np.ndarray, dt: float, omega: float) -> np.ndarray:
    x = np.asarray(x, dtype=float).reshape(4)
    return np.array(
        [
            x[0] + x[2] * dt,
            x[1] + x[3] * dt,
            x[2] - omega * x[3] * dt,
            x[3] + omega * x[2] * dt,
        ]
    )


def motion_jacobian(_x: np.ndarray, dt: float, omega: float) -> np.ndarray:
    return np.array(
        [[1, 0, dt, 0], [0, 1, 0, dt], [0, 0, 1, -omega * dt], [0, 0, omega * dt, 1]],
        dtype=float,
    )


def range_measurement(x: np.ndarray, emitter: np.ndarray) -> float:
    return float(np.linalg.norm(np.asarray(x)[:2] - np.asarray(emitter)[:2]))


def range_jacobian(x: np.ndarray, emitter: np.ndarray) -> np.ndarray:
    delta = np.asarray(x)[:2] - np.asarray(emitter)[:2]
    distance = max(float(np.linalg.norm(delta)), np.finfo(float).eps)
    return np.array([delta[0] / distance, delta[1] / distance, 0.0, 0.0])


class SimulationModel:
    def __init__(self, cfg, data):
        self.cfg, self.data = cfg, data

    def propagate(self, x, index):
        return (
            motion(x, self.cfg.dt, self.cfg.omega),
            motion_jacobian(x, self.cfg.dt, self.cfg.omega),
            self.cfg.q,
        )

    def linearize(self, x, index):
        emitters = self.data["emitter_positions"]
        residual = self.data["toa_measurements"][:, index] - np.array(
            [range_measurement(x, e) for e in emitters.T]
        )
        H = np.vstack([range_jacobian(x, e) for e in emitters.T])
        return residual, H, np.eye(len(residual)) * self.cfg.r
