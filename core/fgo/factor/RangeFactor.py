import numpy as np
from .factor import Factor, whiten, robust_sqrt_weight
from config.config import range_measurement, range_jacobian


class RangeFactor(Factor):
	def evaluate(self):
		measurement = self.z["range"]
		emitter = self.z["emitter"]
		residual = measurement - range_measurement(self.states[0].value, emitter)
		jacobian = -range_jacobian(self.states[0].value, emitter)
		weight = robust_sqrt_weight(
			abs(residual) * np.sqrt(float(self.omega[0, 0])),
			self.z["loss_type"], self.z["loss_delta"])
		self.A, self.b = whiten(weight * jacobian[None, :], weight * residual, self.omega)
		return self

__all__ = ["RangeFactor"]
