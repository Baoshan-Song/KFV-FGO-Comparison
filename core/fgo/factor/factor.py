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
	value = abs(float(residual))
	kind = loss_type.lower()
	if kind == "huber": return 1.0 if value <= delta else np.sqrt(delta / value)
	if kind == "cauchy": return np.sqrt(1.0 / (1.0 + (value / delta) ** 2))
	if kind == "tukey": return np.sqrt((1.0 - (value / delta) ** 2) ** 2) if value <= delta else 0.0
	if kind == "none": return 1.0
	raise ValueError(f"Unknown loss type: {loss_type}")

__all__ = ["Factor"]
