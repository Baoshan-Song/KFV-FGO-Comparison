import numpy as np


def _predict_measurements(x, measurements, emitters, h, H):
	values = np.array([h(x, emitter) for emitter in emitters.T])
	jacobian = np.vstack([H(x, emitter) for emitter in emitters.T])
	return measurements - values, jacobian


def ekf(x_prev, p_prev, dt, omega, f, F, Q, measurements, emitters, h, H, r_base):
    from .miekf import iterated_update
    return iterated_update(x_prev, p_prev, dt, omega, f, F, Q, measurements,
                           emitters, h, H, r_base, 1, 0.)

__all__ = ["ekf"]
