"""Fixed-scale measurement losses and comparable per-iteration diagnostics."""
import numpy as np


def robust_weights(residual, kind='none', delta=2.):
    u = np.asarray(residual, dtype=float)
    if not np.isfinite(delta) or delta <= 0:
        raise ValueError('robust_delta must be finite and positive')
    a = np.abs(u)
    kind = kind.lower()
    if kind == 'none': return np.ones_like(u)
    if kind == 'huber': return np.minimum(1., delta / np.maximum(a, np.finfo(float).tiny))
    if kind == 'cauchy': return (delta / np.hypot(delta, u))**2
    if kind == 'tukey': return (1. - np.minimum(a / delta, 1.)**2)**2
    raise ValueError(f'Unknown loss type: {kind}')


def loss_cost(residual, kind='none', delta=2.):
    """Sum rho(u), with rho(u)=u^2/2 for ordinary least squares."""
    u = np.asarray(residual, dtype=float)
    if not np.isfinite(delta) or delta <= 0:
        raise ValueError('robust_delta must be finite and positive')
    a = np.abs(u)
    kind = kind.lower()
    if kind == 'none': v = .5*u*u
    elif kind == 'huber': v = np.where(a <= delta, .5*u*u, delta*(a-.5*delta))
    elif kind == 'cauchy': v = delta**2*(np.log(np.hypot(delta, u))-np.log(delta))
    elif kind == 'tukey': v = delta**2/6*(1-(1-np.minimum(a/delta, 1.)**2)**3)
    else: raise ValueError(f'Unknown loss type: {kind}')
    return float(np.sum(v))


def continuation_start(residual, kind, target, enabled):
    """Explicit optional broad-to-target scale for redescending losses.

    Cauchy starts at twice the median absolute whitened residual. Tukey
    starts at four times the maximum, keeping its first-step weights at
    least (15/16)^2. These are explicit initialization choices, not
    statistical gates. The nominal measurement covariance never changes.
    Fixed-scale experiments are retained with enabled=False.
    """
    if not enabled or kind not in ('cauchy', 'tukey'):
        return float(target)
    values = np.abs(np.asarray(residual, dtype=float))
    values = values[np.isfinite(values)]
    if not values.size: return float(target)
    broad = 4.*float(np.max(values)) if kind == 'tukey' else 2.*float(np.median(values))
    return max(float(target), broad)


def continuation_scale(start, target, iteration, budget):
    if budget <= 1:
        return start
    if iteration >= budget:
        return target
    return float(target * (start / target)**(1.-(iteration-1)/(budget-1)))
