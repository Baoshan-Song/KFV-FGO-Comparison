"""Four-state simulation propagation and range observations."""

import numpy as np


from config.config import motion, motion_jacobian, range_measurement, range_jacobian


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
