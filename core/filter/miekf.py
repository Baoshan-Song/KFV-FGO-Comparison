"""Fixed-prediction IEKF, arXiv:2511.00306v2, Eqs. (18)--(20)."""
import numpy as np
from .ekf import _predict_measurements
from ..objective import robust_weights, loss_cost


def iterated_update(x_prev, p_prev, dt, omega, f, F, Q, measurements,
                    emitters, h, H, r_base, max_iter, threshold,
                    loss_type='none', delta=2.):
    if max_iter < 1 or r_base <= 0:
        raise ValueError('At least one iteration and positive measurement variance are required')
    x_pred = np.asarray(f(x_prev, dt, omega), dtype=float).reshape(-1)
    jac_f = F(x_prev, dt, omega)
    p_pred = jac_f @ p_prev @ jac_f.T + Q
    x = x_pred.copy()
    identity = np.eye(len(x))
    covariance = np.eye(len(measurements))*float(r_base)
    prior_root = np.linalg.cholesky(p_pred)
    trace = []
    def record(value, iteration):
        raw, _ = _predict_measurements(value, measurements, emitters, h, H)
        prior = np.linalg.solve(prior_root, value-x_pred)
        trace.append(dict(iteration=iteration, state=value.copy(),
            measurement_residual_norm_m=float(np.linalg.norm(raw)),
            objective=.5*float(prior@prior)+loss_cost(raw/np.sqrt(r_base), loss_type, delta)))
    record(x, 0)
    converged = False
    for iteration in range(1, max_iter+1):
        linearization = x.copy()
        residual, jacobian = _predict_measurements(x, measurements, emitters, h, H)
        weights = robust_weights(residual/np.sqrt(r_base), loss_type, delta)
        root = np.sqrt(weights)
        weighted_h = root[:, None]*jacobian
        innovation = weighted_h @ p_pred @ weighted_h.T + covariance
        gain = np.linalg.solve(innovation, weighted_h @ p_pred).T
        # The prediction is fixed throughout this update. This correction is
        # essential even for a linear observation iterated more than once.
        x_next = x_pred + gain @ (root*(residual + jacobian@(x-x_pred)))
        correction = identity - gain @ weighted_h
        p_est = correction @ p_pred @ correction.T + gain @ covariance @ gain.T
        p_est = .5*(p_est+p_est.T)
        step = x_next-x
        x = x_next
        record(x, iteration)
        converged = np.linalg.norm(step)/len(x) < threshold
        if converged: break
    postfit, _ = _predict_measurements(x, measurements, emitters, h, H)
    prefit, _ = _predict_measurements(x_pred, measurements, emitters, h, H)
    return x, p_est, x_pred, p_pred, dict(
        jacobian_all=jacobian, kalman_gain=gain, innovation_covariance=innovation,
        residual=postfit, prefit_residual=prefit, postfit_residual=postfit,
        residual_norm_all=[t['measurement_residual_norm_m'] for t in trace],
        iteration_trace=trace, converged=converged, iterations=iteration,
        linearization_point=linearization, measurement_weights=weights)


def miekf(x_prev, p_prev, dt, omega, f, F, Q, measurements, emitters, h, H,
          r_base, max_iter, threshold):
    return iterated_update(x_prev, p_prev, dt, omega, f, F, Q, measurements,
                           emitters, h, H, r_base, max_iter, threshold)

__all__ = ['miekf']
