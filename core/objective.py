"""Fixed-scale measurement losses and comparable per-iteration diagnostics."""
import numpy as np


def robust_weights(residual, kind='none', delta=2.):
    u = np.asarray(residual, dtype=float)
    if delta <= 0:
        raise ValueError('robust_delta must be positive')
    a = np.abs(u)
    kind = kind.lower()
    if kind == 'none': return np.ones_like(u)
    if kind == 'huber': return np.minimum(1., delta / np.maximum(a, np.finfo(float).tiny))
    if kind == 'cauchy': return 1. / (1. + (u / delta)**2)
    if kind == 'tukey': return np.where(a < delta, (1. - (u / delta)**2)**2, 0.)
    raise ValueError(f'Unknown loss type: {kind}')


def loss_cost(residual, kind='none', delta=2.):
    """Sum rho(u), with rho(u)=u^2/2 for ordinary least squares."""
    u = np.asarray(residual, dtype=float)
    a = np.abs(u)
    kind = kind.lower()
    if kind == 'none': v = .5*u*u
    elif kind == 'huber': v = np.where(a <= delta, .5*u*u, delta*(a-.5*delta))
    elif kind == 'cauchy': v = .5*delta**2*np.log1p((u/delta)**2)
    elif kind == 'tukey': v = np.where(a < delta, delta**2/6*(1-(1-(u/delta)**2)**3), delta**2/6)
    else: raise ValueError(f'Unknown loss type: {kind}')
    return float(np.sum(v))
