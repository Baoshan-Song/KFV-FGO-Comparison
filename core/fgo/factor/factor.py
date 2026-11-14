from dataclasses import dataclass
import numpy as np


@dataclass
class State:
	gid: int
	lid: int
	value: np.ndarray
	status: str = "Forward"


class Factor:
	def __init__(self, states, z=None, omega=None):
		self.states = list(states)
		self.z = z
		self.omega = omega
		self.A = None
		self.b = None
		self.status = ""

	def evaluate(self):
		raise NotImplementedError


def whiten(A, b, omega):
	lower = np.linalg.cholesky(np.asarray(omega, dtype=float))
	return lower.T @ A, lower.T @ np.asarray(b).reshape(-1)


def robust_sqrt_weight(residual, loss_type, delta):
	from core.objective import robust_weights
	return float(np.sqrt(robust_weights(residual, loss_type, delta)))

__all__ = ["Factor"]
