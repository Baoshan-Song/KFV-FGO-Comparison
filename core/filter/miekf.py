import numpy as np
from .ekf import _predict_measurements


def miekf(x_prev, p_prev, dt, omega, f, F, Q, measurements, emitters, h, H,
		  r_base, max_iter, threshold):
	x_pred = f(x_prev, dt, omega)
	p_pred = F(x_prev, dt, omega) @ p_prev @ F(x_prev, dt, omega).T + Q
	x_last = x_pred.copy()
	identity = np.eye(len(x_pred))
	covariance = np.eye(len(measurements)) * r_base
	for _ in range(max_iter):
		residual, jacobian = _predict_measurements(x_last, measurements, emitters, h, H)
		innovation = jacobian @ p_pred @ jacobian.T + covariance
		gain = np.linalg.solve(innovation.T, (p_pred @ jacobian.T).T).T
		x_est = x_last + gain @ residual
		p_est = (identity - gain @ jacobian) @ p_pred
		if np.linalg.norm(x_est - x_last) / len(x_est) < threshold:
			break
		x_last = x_est
	return x_est, p_est, x_pred, p_pred, {
		"jacobian_all": jacobian, "residual_norm_all": [],
		"kalman_gain": gain, "innovation_covariance": innovation,
		"residual": residual,
	}

__all__ = ["miekf"]
