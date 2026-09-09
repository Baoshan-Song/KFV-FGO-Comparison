import numpy as np
from .factor import Factor, whiten


class AutoDiffFactor(Factor):
	def __init__(self, states, z, omega, error_func, delta=1e-4):
		super().__init__(states, z, omega)
		self.error_func = error_func
		self.delta = delta

	def evaluate(self):
		x = self.states[0].value
		residual = np.asarray(self.error_func(x, self.z), dtype=float).reshape(-1)
		jacobian = np.zeros((len(residual), len(x)))
		for index in range(len(x)):
			plus, minus = x.copy(), x.copy()
			plus[index] += self.delta
			minus[index] -= self.delta
			jacobian[:, index] = (self.error_func(plus, self.z) - self.error_func(minus, self.z)) / (2 * self.delta)
		self.A, self.b = whiten(jacobian, residual, self.omega)
		return self

__all__ = ["AutoDiffFactor"]
