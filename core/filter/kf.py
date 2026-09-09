import numpy as np


def kf(x_prev, p_prev, F, Q, B, U, measurement, H, R):
	x_pred = F @ x_prev + B @ U
	p_pred = F @ p_prev @ F.T + Q
	gain = np.linalg.solve(H @ p_pred @ H.T + R, H @ p_pred).T
	x_est = x_pred + gain @ (measurement - H @ x_pred)
	p_est = (np.eye(len(x_pred)) - gain @ H) @ p_pred
	return x_est, p_est, x_pred, p_pred

__all__ = ["kf"]
