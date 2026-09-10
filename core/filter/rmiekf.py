import numpy as np
from .ekf import _predict_measurements


def compute_weight_matlab(r, loss_type, delta):
    """严格对齐 MATLAB compute_weight 子函数"""
    abs_r = np.abs(r)
    loss = str(loss_type).lower()

    if loss == 'huber':
        if abs_r <= delta:
            return 1.0
        else:
            return float(delta / abs_r)
    elif loss == 'cauchy':
        return float(1.0 / (1.0 + (r / delta) ** 2))
    elif loss == 'tukey':
        if abs_r <= delta:
            return float((1.0 - (r / delta) ** 2) ** 2)
        else:
            return 0.0
    elif loss == 'none':
        return 1.0
    else:
        raise ValueError(f"Unknown loss type: {loss_type}")


def rmiekf(x_prev, p_prev, dt, omega, f, F, Q, measurements, emitters, h, H,
           r_base, max_iter, threshold, loss_type, delta):
    
    # 强制将一维状态和测量转换为列向量/一维标准格式，防止广播 Bug
    x_prev = np.asarray(x_prev).flatten()
    measurements = np.asarray(measurements).flatten()

    # 1. Prediction
    x_pred = f(x_prev, dt, omega).flatten()
    p_pred = F(x_prev, dt, omega) @ p_prev @ F(x_prev, dt, omega).T + Q

    x_est = x_pred.copy()
    p_est = p_pred.copy()
    identity = np.eye(len(x_est))
    x_last = x_est.copy()

    num_emitters = len(measurements)
    # MATLAB: R = R_base * eye(num_emitters);
    R = np.eye(num_emitters) * float(r_base)

    gain = None
    innovation = None
    jacobian = None
    residual = None

    # 2. Iterative Measurement Update
    for _ in range(max_iter):
        # 确保残差和雅可比计算与 MATLAB 格式一致
        residual, jacobian = _predict_measurements(x_last, measurements, emitters, h, H)
        residual = residual.flatten()  # 确保是 1D 向量 (m,)

        # MATLAB: r_vec = abs(y) ./ sqrt(diag(R));
        r_diag = np.sqrt(np.diag(R))
        r_vec = np.abs(residual) / r_diag

        # MATLAB: w_vec = arrayfun(@(r) compute_weight(r, loss_type, delta), r_vec);
        w_vec = np.array([compute_weight_matlab(r, loss_type, delta) for r in r_vec])

        # MATLAB: W = diag(w_vec); H_r = W^(1/2) * H_all; y_r = W^(1/2) * y;
        root_w = np.sqrt(w_vec)
        H_r = root_w[:, None] * jacobian  # m x n
        y_r = root_w * residual            # m x 1

        # MATLAB: R_r = R;
        R_r = R.copy()

        # MATLAB: R = diag(1 ./ w_vec) * R;
        # 增加防御：防止 1/0 出现 Inf 导致协方差破坏
        inv_w_vec = np.where(w_vec > 1e-12, 1.0 / w_vec, 1e12)
        R = np.diag(inv_w_vec) @ R

        # MATLAB: S_r = H_r * P_pred * H_r' + R_r;
        innovation = H_r @ p_pred @ H_r.T + R_r

        # MATLAB: K = P_pred * H_r' / S_r;
        # 严格等价于 MATLAB 的右除 / 运算
        gain = np.linalg.solve(innovation, H_r @ p_pred).T

        # MATLAB: X_est = X_last_est + K * y_r;
        x_est = x_last + gain @ y_r

        # MATLAB: P_est = (I - K * H_r) * P_pred;
        p_est = (identity - gain @ H_r) @ p_pred

        # MATLAB: if norm(X_est-X_last_est)/length(X_est) < thres
        if np.linalg.norm(x_est - x_last) / len(x_est) < threshold:
            break

        # MATLAB: X_last_est = X_est; (使用 .copy() 保证独立内存)
        x_last = x_est.copy()

    debug_info = {
        "jacobian_all": jacobian,
        "residual_norm_all": [],
        "kalman_gain": gain,
        "innovation_covariance": innovation,
        "residual": residual,
    }

    return x_est, p_est, x_pred, p_pred, debug_info


__all__ = ["rmiekf"]


# import numpy as np
# from .ekf import _predict_measurements
# from .rekf import _weight


# def rmiekf(x_prev, p_prev, dt, omega, f, F, Q, measurements, emitters, h, H,
# 		   r_base, max_iter, threshold, loss_type, delta):
# 	x_pred = f(x_prev, dt, omega)
# 	p_pred = F(x_prev, dt, omega) @ p_prev @ F(x_prev, dt, omega).T + Q
# 	x_last = x_pred.copy()
# 	identity = np.eye(len(x_pred))
# 	# MATLAB keeps and updates R across iterations.
# 	covariance = np.eye(len(measurements)) * r_base
# 	for _ in range(max_iter):
# 		residual, jacobian = _predict_measurements(x_last, measurements, emitters, h, H)
# 		normalized = np.abs(residual) / np.sqrt(np.diag(covariance))
# 		weights = np.array([_weight(value, loss_type, delta) for value in normalized])
# 		root = np.sqrt(weights)
# 		weighted_jacobian = root[:, None] * jacobian
# 		weighted_residual = root * residual
# 		robust_covariance = covariance.copy()
# 		innovation = weighted_jacobian @ p_pred @ weighted_jacobian.T + robust_covariance
# 		gain = np.linalg.solve(innovation.T, (p_pred @ weighted_jacobian.T).T).T
# 		x_est = x_last + gain @ weighted_residual
# 		p_est = (identity - gain @ weighted_jacobian) @ p_pred
# 		covariance = np.diag(1.0 / weights) @ covariance
# 		if np.linalg.norm(x_est - x_last) / len(x_est) < threshold:
# 			break
# 		x_last = x_est
# 	return x_est, p_est, x_pred, p_pred, {
# 		"jacobian_all": jacobian, "residual_norm_all": [],
# 		"kalman_gain": gain, "innovation_covariance": innovation,
# 		"residual": residual,
# 	}

# __all__ = ["rmiekf"]
