"""MATLAB-compatible KFV recurrences, including its modified iterated update."""

from time import perf_counter

import numpy as np


def _weight(residual, loss_type, delta):
    value = abs(float(residual))
    kind = loss_type.lower()
    if kind == "none":
        return 1.0
    if kind == "huber":
        return 1.0 if value <= delta else delta / value
    if kind == "cauchy":
        return 1.0 / (1 + (value / delta) ** 2)
    if kind == "tukey":
        return (1 - (value / delta) ** 2) ** 2 if value <= delta else 0.0
    raise ValueError(f"Unknown robust kernel: {loss_type}")


def measurement_covariance(value, n):
    array = np.asarray(value, dtype=float)
    if array.ndim == 0:
        return np.eye(n) * array
    if array.ndim == 1:
        return np.diag(array)
    if array.shape != (n, n):
        raise ValueError("Measurement covariance dimension mismatch")
    return array.copy()


def _linearized_measurements(x, z, emitters, h, H):
    return (
        np.asarray(z) - np.array([h(x, e) for e in emitters.T]),
        np.vstack([H(x, e) for e in emitters.T]),
    )


def filter_step(
    x,
    p,
    predicted,
    transition,
    q,
    observe,
    *,
    mode="EKF",
    max_iter=5,
    threshold=1e-6,
    loss_type="huber",
    delta=1,
    dynamic_observation=True,
    trace=False,
):
    p_pred = transition @ p @ transition.T + q
    last = predicted.copy()
    robust = mode in ("rEKF", "riEKF")
    limit = max_iter if mode in ("iEKF", "riEKF") else 1
    if mode not in ("EKF", "iEKF", "rEKF", "riEKF"):
        raise ValueError(f"Unsupported KFV mode: {mode}")
    iterations = []
    covariance = None
    started = perf_counter()
    for iteration in range(limit):
        residual, jacobian, evaluated = observe(last)[:3]
        residual = np.asarray(residual).reshape(-1)
        if iteration == 0 or dynamic_observation:
            covariance = measurement_covariance(evaluated, len(residual))
        weights = np.ones(len(residual))
        if robust:
            standardized = residual / np.sqrt(np.diag(covariance))
            weights = np.array([_weight(v, loss_type, delta) for v in standardized])
        root = np.sqrt(weights)
        h = root[:, None] * jacobian
        y = root * residual
        innovation = h @ p_pred @ h.T + covariance
        gain = np.linalg.solve(innovation.T, (p_pred @ h.T).T).T
        estimate = last + gain @ y
        posterior = (np.eye(len(x)) - gain @ h) @ p_pred
        if trace:
            iterations.append(
                {
                    "linearization_state": last.copy(),
                    "residual": residual.copy(),
                    "H": jacobian.copy(),
                    "R": covariance.copy(),
                    "weights": weights.copy(),
                    "X": estimate.copy(),
                    "P": posterior.copy(),
                }
            )
        # Legacy simulation riEKF accumulates R inflation; real callback resets R.
        if mode == "riEKF" and not dynamic_observation:
            with np.errstate(divide="ignore", invalid="ignore"):
                covariance = np.diag(np.diag(covariance) / weights)
        difference = np.linalg.norm(estimate - last) / len(x)
        last = estimate
        if difference < threshold:
            break
    if not np.all(np.isfinite(estimate)) or not np.all(np.isfinite(posterior)):
        raise FloatingPointError("KFV update produced non-finite values")
    # MATLAB rEKF exposes the unweighted S in debug, although it solves with S_r.
    reported_innovation = (
        jacobian @ p_pred @ jacobian.T + covariance if mode == "rEKF" else innovation
    )
    debug = {
        "jacobian_all": jacobian,
        "residual": residual,
        "kalman_gain": gain,
        "innovation_covariance": reported_innovation,
        "weighted_innovation_covariance": innovation,
        "iterations": iterations,
        "residual_norm_all": residual**2,
        "iteration_count": iteration + 1,
        "update_time": perf_counter() - started,
        "X_pred": predicted.copy(),
        "P_pred": p_pred,
        "P": posterior.copy(),
    }
    return estimate, posterior, predicted, p_pred, debug


def _legacy(
    x,
    p,
    dt,
    omega,
    f,
    F,
    Q,
    z,
    emitters,
    h,
    H,
    R,
    mode,
    max_iter=1,
    threshold=1e-6,
    loss_type="none",
    delta=1,
    observation_function=None,
):
    def observe(state):
        if observation_function is not None:
            return observation_function(state, z)
        residual, jacobian = _linearized_measurements(state, z, emitters, h, H)
        return residual, jacobian, R

    return filter_step(
        x,
        p,
        f(x, dt, omega),
        F(x, dt, omega),
        Q,
        observe,
        mode=mode,
        max_iter=max_iter,
        threshold=threshold,
        loss_type=loss_type,
        delta=delta,
        dynamic_observation=observation_function is not None,
    )


def ekf(x, p, dt, omega, f, F, Q, z, emitters, h, H, R, observation_function=None):
    return _legacy(
        x,
        p,
        dt,
        omega,
        f,
        F,
        Q,
        z,
        emitters,
        h,
        H,
        R,
        "EKF",
        observation_function=observation_function,
    )


def miekf(
    x,
    p,
    dt,
    omega,
    f,
    F,
    Q,
    z,
    emitters,
    h,
    H,
    R,
    max_iter,
    threshold,
    observation_function=None,
):
    return _legacy(
        x,
        p,
        dt,
        omega,
        f,
        F,
        Q,
        z,
        emitters,
        h,
        H,
        R,
        "iEKF",
        max_iter,
        threshold,
        observation_function=observation_function,
    )


def rekf(
    x,
    p,
    dt,
    omega,
    f,
    F,
    Q,
    z,
    emitters,
    h,
    H,
    R,
    loss_type,
    delta,
    observation_function=None,
):
    return _legacy(
        x,
        p,
        dt,
        omega,
        f,
        F,
        Q,
        z,
        emitters,
        h,
        H,
        R,
        "rEKF",
        loss_type=loss_type,
        delta=delta,
        observation_function=observation_function,
    )


def rmiekf(
    x,
    p,
    dt,
    omega,
    f,
    F,
    Q,
    z,
    emitters,
    h,
    H,
    R,
    max_iter,
    threshold,
    loss_type,
    delta,
    observation_function=None,
):
    return _legacy(
        x,
        p,
        dt,
        omega,
        f,
        F,
        Q,
        z,
        emitters,
        h,
        H,
        R,
        "riEKF",
        max_iter,
        threshold,
        loss_type,
        delta,
        observation_function,
    )


def kf(x_prev, p_prev, F, Q, B, U, measurement, H, R):
    x_pred = F @ x_prev + B @ U
    p_pred = F @ p_prev @ F.T + Q
    gain = np.linalg.solve(H @ p_pred @ H.T + R, H @ p_pred).T
    return (
        x_pred + gain @ (measurement - H @ x_pred),
        (np.eye(len(x_pred)) - gain @ H) @ p_pred,
        x_pred,
        p_pred,
    )
