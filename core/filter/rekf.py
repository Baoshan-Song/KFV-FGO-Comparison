"""Single-linearization robust update using fixed measurement noise."""
from .miekf import iterated_update
from ..objective import robust_weights


def _weight(residual, loss_type, delta):
    return float(robust_weights(residual, loss_type, delta))


def rekf(x_prev, p_prev, dt, omega, f, F, Q, measurements, emitters, h, H,
         r_base, loss_type, delta):
    return iterated_update(x_prev, p_prev, dt, omega, f, F, Q, measurements,
                           emitters, h, H, r_base, 1, 0., loss_type, delta)

__all__ = ['rekf']
