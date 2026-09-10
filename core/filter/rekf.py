import numpy as np
from .ekf import _predict_measurements


def _weight(residual, loss_type, delta):
	value = abs(float(residual))
	kind = loss_type.lower()
	if kind == "huber": return 1.0 if value <= delta else delta / value
	if kind == "cauchy": return 1.0 / (1.0 + (value / delta) ** 2)
	if kind == "tukey": return (1.0 - (value / delta) ** 2) ** 2 if value <= delta else 0.0
	if kind == "none": return 1.0
	raise ValueError(f"Unknown loss type: {loss_type}")


def rekf(x_prev, p_prev, dt, omega, f, F, Q, measurements, emitters, h, H,
		 r_base, loss_type, delta):
	x_pred = f(x_prev, dt, omega)
	p_pred = F(x_prev, dt, omega) @ p_prev @ F(x_prev, dt, omega).T + Q
	residual, jacobian = _predict_measurements(x_pred, measurements, emitters, h, H)
	weights = np.array([_weight(value / np.sqrt(r_base), loss_type, delta) for value in residual])
	root = np.sqrt(weights)
	weighted_jacobian = root[:, None] * jacobian
	weighted_residual = root * residual
	covariance = np.eye(len(measurements)) * r_base
	raw_innovation = jacobian @ p_pred @ jacobian.T + covariance
	# innovation = weighted_jacobian @ p_pred @ weighted_jacobian.T + covariance
	# gain = np.linalg.solve(innovation.T, (p_pred @ weighted_jacobian.T).T).T
	robust_innovation = weighted_jacobian @ p_pred @ weighted_jacobian.T + covariance
	gain = np.linalg.solve(robust_innovation, weighted_jacobian @ p_pred).T
	x_est = x_pred + gain @ weighted_residual
	p_est = (np.eye(len(x_pred)) - gain @ weighted_jacobian) @ p_pred
	return x_est, p_est, x_pred, p_pred, {
		"jacobian_all": jacobian, "residual_norm_all": residual ** 2,
		"kalman_gain": gain, "innovation_covariance": raw_innovation,
		"residual": residual,
	}

__all__ = ["rekf"]
