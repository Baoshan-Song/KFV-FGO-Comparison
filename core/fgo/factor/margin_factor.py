import numpy as np
from .factor import Factor


class MarginFactor(Factor):
    def __init__(self, states, A, b, x0):
        super().__init__(states, None, None)
        self.A0 = np.asarray(A, dtype=float).copy()
        self.b0 = np.asarray(b, dtype=float).reshape(-1).copy()
        self.x0 = np.asarray(x0, dtype=float).reshape(-1).copy()

    def evaluate(self):
        x = np.concatenate([
            state.value for state in self.states
        ])

        dx = x - self.x0

        self.A = self.A0
        self.b = self.b0 - self.A0 @ dx

        return self

__all__ = ["MarginFactor"]
