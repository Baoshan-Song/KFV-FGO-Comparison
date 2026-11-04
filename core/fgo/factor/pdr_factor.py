import numpy as np
from .factor import Factor, whiten


class PdrFactor(Factor):
	def evaluate(self):
		residual = self.z - (self.states[1].value - self.states[0].value)
		self.A, self.b = whiten(np.hstack((-np.eye(len(self.states[0].value)), np.eye(len(self.states[1].value)))), residual, self.omega)
		return self

__all__ = ["PdrFactor"]
