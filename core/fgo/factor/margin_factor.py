import numpy as np
from .factor import Factor


class MarginFactor(Factor):
    def __init__(self, states, A, b, x0, *, matlab_compat=False):
        super().__init__(states, None, None)
        self.A0 = A
        self.b0 = np.asarray(b).reshape(-1)
        self.x0 = np.asarray(x0).reshape(-1)
        self.matlab_compat = matlab_compat

    def evaluate(self):
        if self.matlab_compat:
            # MATLAB margin_factor.evaluate is a no-op. This intentionally
            # reproduces its recurrence, not a fixed-mean MAP prior.
            self.A, self.b = self.A0, self.b0
            return self

        x = np.concatenate([
            state.value for state in self.states
        ])

        dx = x - self.x0

        self.A = self.A0
        self.b = self.b0 - self.A0 @ dx

        return self

__all__ = ["MarginFactor"]
