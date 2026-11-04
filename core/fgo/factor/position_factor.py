import numpy as np
from .factor import Factor, whiten


class PositionFactor(Factor):
	def evaluate(self):
		# The graph solves A * delta = b and then performs x += delta.
		self.A, self.b = whiten(np.eye(len(self.states[0].value)), self.z - self.states[0].value, self.omega)
		return self

__all__ = ["PositionFactor"]
