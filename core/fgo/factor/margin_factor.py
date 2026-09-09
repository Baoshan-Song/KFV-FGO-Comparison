import numpy as np
from .factor import Factor


class MarginFactor(Factor):
	def __init__(self, states, A, b, omega=None):
		super().__init__(states, None, omega)
		self.A = A
		self.b = np.asarray(b).reshape(-1)

	def evaluate(self):
		return self

__all__ = ["MarginFactor"]
