import numpy as np
from .factor import Factor, whiten
from config.config import motion, motion_jacobian


class PropagateFactor(Factor):
	def __init__(self, states, config):
		super().__init__(states)
		self.config = config

	def evaluate(self):
		cfg = self.config
		jacobian = motion_jacobian(self.states[0].value, cfg.dt, cfg.omega)
		residual = self.states[1].value - motion(self.states[0].value, cfg.dt, cfg.omega)
		self.A, self.b = whiten(np.hstack((jacobian, -np.eye(4))), residual, np.linalg.inv(cfg.q))
		return self

__all__ = ["PropagateFactor"]
