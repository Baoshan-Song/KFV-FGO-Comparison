import numpy as np
from .AutoDiffFactor import AutoDiffFactor


class AutoDiffRobustFactor(AutoDiffFactor):
	def __init__(self, states, z, error_func, kernel_func):
		super().__init__(states, z, np.eye(1), error_func)
		self.kernel_func = kernel_func

	def evaluate(self):
		x = self.states[0].value
		residual = np.asarray(self.kernel_func(self.error_func(x, self.z))).reshape(-1)
		jacobian = np.zeros((len(residual), len(x)))
		for index in range(len(x)):
			plus, minus = x.copy(), x.copy()
			plus[index] += self.delta
			minus[index] -= self.delta
			jacobian[:, index] = (self.kernel_func(self.error_func(plus, self.z)) -
								  self.kernel_func(self.error_func(minus, self.z))) / (2 * self.delta)
		self.A, self.b = jacobian, residual
		return self

__all__ = ["AutoDiffRobustFactor"]
