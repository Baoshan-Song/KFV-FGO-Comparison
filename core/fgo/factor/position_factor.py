import numpy as np
from .factor import Factor, whiten


class PositionFactor(Factor):
	def evaluate(self):
		self.A, self.b = whiten(-np.eye(4), self.z - self.states[0].value, self.omega)
		return self

__all__ = ["PositionFactor"]
