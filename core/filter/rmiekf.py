"""Robust IEKF with a Gaussian prior and fixed-scale measurement IRLS."""
from .miekf import iterated_update


def rmiekf(x_prev, p_prev, dt, omega, f, F, Q, measurements, emitters, h, H,
           r_base, max_iter, threshold, loss_type, delta):
    return iterated_update(x_prev, p_prev, dt, omega, f, F, Q, measurements,
                           emitters, h, H, r_base, max_iter, threshold, loss_type, delta)

__all__ = ['rmiekf']
