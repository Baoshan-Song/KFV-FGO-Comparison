from time import perf_counter
import numpy as np


def _predict_measurements(x, measurements, emitters, h, H):
	values = np.array([h(x, emitter) for emitter in emitters.T])
	jacobian = np.vstack([H(x, emitter) for emitter in emitters.T])
	return measurements - values, jacobian


def ekf(x_prev, p_prev, dt, omega, f, F, Q, measurements, emitters, h, H, r_base):
	prediction_started = perf_counter()
	x_pred = f(x_prev, dt, omega)
	p_pred = F(x_prev, dt, omega) @ p_prev @ F(x_prev, dt, omega).T + Q
	prediction_time = perf_counter() - prediction_started

	residual, jacobian = _predict_measurements(x_pred, measurements, emitters, h, H)
	covariance = np.eye(len(measurements)) * r_base
	innovation = jacobian @ p_pred @ jacobian.T + covariance
	gain = np.linalg.solve(innovation.T, (p_pred @ jacobian.T).T).T
	x_est = x_pred + gain @ residual
	p_est = (np.eye(len(x_pred)) - gain @ jacobian) @ p_pred
	debug_info = {
		"prediction_time": prediction_time,
		"jacobian_all": jacobian,
		"residual_norm_all": residual ** 2,
		"kalman_gain": gain,
		"innovation_covariance": innovation,
		"residual": residual,
	}
	return x_est, p_est, x_pred, p_pred, debug_info

__all__ = ["ekf"]
