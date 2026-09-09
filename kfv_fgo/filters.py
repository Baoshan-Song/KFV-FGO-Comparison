from __future__ import annotations

from time import perf_counter
import numpy as np


def _weight(residual, loss_type, delta):
    r = abs(float(residual))
    kind = loss_type.lower()
    if kind == "huber": return 1.0 if r <= delta else delta / r
    if kind == "cauchy": return 1.0 / (1.0 + (r / delta) ** 2)
    if kind == "tukey": return (1.0 - (r / delta) ** 2) ** 2 if r <= delta else 0.0
    if kind == "none": return 1.0
    raise ValueError(f"Unknown loss type: {loss_type}")


def _predict(x, p, dt, omega, f, F, Q):
    started = perf_counter()
    x_pred = f(x, dt, omega)
    p_pred = F(x, dt, omega) @ p @ F(x, dt, omega).T + Q
    return x_pred, p_pred, perf_counter() - started


def _linearized_measurements(x, z, emitters, h, H):
    values = np.array([h(x, emitter) for emitter in emitters.T])
    jacobian = np.vstack([H(x, emitter) for emitter in emitters.T])
    return z - values, jacobian


def ekf(x_prev, p_prev, dt, omega, f, F, Q, measurements, emitters, h, H, r_base):
    x_pred, p_pred, prediction_time = _predict(x_prev, p_prev, dt, omega, f, F, Q)
    residual, jacobian = _linearized_measurements(x_pred, measurements, emitters, h, H)
    covariance = np.eye(len(measurements)) * r_base
    innovation = jacobian @ p_pred @ jacobian.T + covariance
    gain = np.linalg.solve(innovation.T, (p_pred @ jacobian.T).T).T
    x_est = x_pred + gain @ residual
    p_est = (np.eye(len(x_pred)) - gain @ jacobian) @ p_pred
    return x_est, p_est, x_pred, p_pred, {"prediction_time": prediction_time, "jacobian_all": jacobian,
        "residual_norm_all": residual ** 2, "kalman_gain": gain, "innovation_covariance": innovation,
        "residual": residual}


def miekf(x_prev, p_prev, dt, omega, f, F, Q, measurements, emitters, h, H, r_base, max_iter, threshold):
    x_pred, p_pred, prediction_time = _predict(x_prev, p_prev, dt, omega, f, F, Q)
    x_last = x_pred.copy()
    identity = np.eye(len(x_pred))
    covariance = np.eye(len(measurements)) * r_base
    for _ in range(max_iter):
        residual, jacobian = _linearized_measurements(x_last, measurements, emitters, h, H)
        innovation = jacobian @ p_pred @ jacobian.T + covariance
        gain = np.linalg.solve(innovation.T, (p_pred @ jacobian.T).T).T
        x_est = x_last + gain @ residual
        p_est = (identity - gain @ jacobian) @ p_pred
        if np.linalg.norm(x_est - x_last) / len(x_est) < threshold:
            break
        x_last = x_est
    return x_est, p_est, x_pred, p_pred, {"prediction_time": prediction_time, "jacobian_all": jacobian,
        "residual_norm_all": [], "kalman_gain": gain, "innovation_covariance": innovation, "residual": residual}


def rekf(x_prev, p_prev, dt, omega, f, F, Q, measurements, emitters, h, H, r_base, loss_type, delta):
    x_pred, p_pred, prediction_time = _predict(x_prev, p_prev, dt, omega, f, F, Q)
    residual, jacobian = _linearized_measurements(x_pred, measurements, emitters, h, H)
    weights = np.array([_weight(value / np.sqrt(r_base), loss_type, delta) for value in residual])
    root = np.sqrt(weights)
    weighted_jacobian, weighted_residual = root[:, None] * jacobian, root * residual
    covariance = np.eye(len(measurements)) * r_base
    innovation = weighted_jacobian @ p_pred @ weighted_jacobian.T + covariance
    gain = np.linalg.solve(innovation.T, (p_pred @ weighted_jacobian.T).T).T
    x_est = x_pred + gain @ weighted_residual
    p_est = (np.eye(len(x_pred)) - gain @ weighted_jacobian) @ p_pred
    return x_est, p_est, x_pred, p_pred, {"prediction_time": prediction_time, "jacobian_all": jacobian,
        "residual_norm_all": residual ** 2, "kalman_gain": gain, "innovation_covariance": innovation,
        "residual": residual}


def rmiekf(x_prev, p_prev, dt, omega, f, F, Q, measurements, emitters, h, H, r_base, max_iter, threshold, loss_type, delta):
    x_pred, p_pred, prediction_time = _predict(x_prev, p_prev, dt, omega, f, F, Q)
    x_last = x_pred.copy()
    identity = np.eye(len(x_pred))
    base_r = np.eye(len(measurements)) * r_base
    for _ in range(max_iter):
        residual, jacobian = _linearized_measurements(x_last, measurements, emitters, h, H)
        weights = np.array([_weight(value / np.sqrt(r_base), loss_type, delta) for value in residual])
        root = np.sqrt(weights)
        weighted_jacobian, weighted_residual = root[:, None] * jacobian, root * residual
        innovation = weighted_jacobian @ p_pred @ weighted_jacobian.T + base_r
        gain = np.linalg.solve(innovation.T, (p_pred @ weighted_jacobian.T).T).T
        x_est = x_last + gain @ weighted_residual
        p_est = (identity - gain @ weighted_jacobian) @ p_pred
        if np.linalg.norm(x_est - x_last) / len(x_est) < threshold:
            break
        x_last = x_est
    return x_est, p_est, x_pred, p_pred, {"prediction_time": prediction_time, "jacobian_all": jacobian,
        "residual_norm_all": [], "kalman_gain": gain, "innovation_covariance": innovation, "residual": residual}


def kf(x_prev, p_prev, F, Q, B, U, measurement, H, R):
    x_pred = F @ x_prev + B @ U
    p_pred = F @ p_prev @ F.T + Q
    gain = np.linalg.solve(H @ p_pred @ H.T + R, H @ p_pred).T
    x_est = x_pred + gain @ (measurement - H @ x_pred)
    return x_est, (np.eye(len(x_pred)) - gain @ H) @ p_pred, x_pred, p_pred
