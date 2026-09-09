from __future__ import annotations

import numpy as np
from .factors import MarginFactor


class FactorGraph:
    def __init__(self, config):
        self.config = config; self.states = []; self.factors = []; self.win_size = 0
        self.J = self.r = None; self.latest_information_matrix = None
        self.ls_time = self.margin_time = self.margin_meas_time = 0.0
        self.residual_norm_all = []

    def add_state(self, state):
        state.lid = sum(s.status != "Margin" for s in self.states) + 1
        self.states.append(state); self.win_size += 1; return self

    def add_factor(self, factor): self.factors.append(factor); return self

    @property
    def active_states(self): return [state for state in self.states if state.status != "Margin"]

    def normal_equation(self):
        active = [factor.evaluate() for factor in self.factors if factor.status != "Margin"]
        n = len(self.active_states); size = len(self.active_states[0].value)
        offsets = {state.gid: index for index, state in enumerate(self.active_states)}
        rows = sum(len(factor.b) for factor in active)
        J = np.zeros((rows, n * size)); r = np.zeros(rows); row = 0
        for factor in active:
            count = len(factor.b); r[row:row + count] = factor.b
            for index, state in enumerate(factor.states):
                if state.status != "Margin":
                    column = offsets[state.gid] * size
                    J[row:row + count, column:column + size] = factor.A[:, index * size:(index + 1) * size]
            row += count
        self.J, self.r = J, r; return self

    def estimate(self):
        for _ in range(self.config.max_iteration):
            self.normal_equation()
            delta, *_ = np.linalg.lstsq(self.J, -self.r, rcond=None)
            self.latest_information_matrix = self.J.T @ self.J
            for index, state in enumerate(self.active_states):
                size = len(state.value); state.value += delta[index * size:(index + 1) * size]
            self.residual_norm_all.append(np.r_[self.active_states[0].value, np.linalg.norm(self.r)])
            if np.linalg.norm(delta) / len(delta) < self.config.threshold_iteration: break
        return self

    def marginalize(self, gids):
        remove = {gids} if np.isscalar(gids) else set(gids)
        self.normal_equation(); active = self.active_states; size = len(active[0].value)
        remove_columns = [i for i, state in enumerate(active) if state.gid in remove]
        keep_columns = [i for i, state in enumerate(active) if state.gid not in remove]
        if remove_columns and keep_columns:
            rm = np.concatenate([np.arange(i * size, (i + 1) * size) for i in remove_columns])
            keep = np.concatenate([np.arange(i * size, (i + 1) * size) for i in keep_columns])
            H = self.J.T @ self.J; b = self.J.T @ self.r
            H11, H12 = H[np.ix_(rm, rm)], H[np.ix_(rm, keep)]
            H22, b1, b2 = H[np.ix_(keep, keep)], b[rm], b[keep]
            Hm = H22 - H[np.ix_(keep, rm)] @ np.linalg.solve(H11, H12)
            bm = b2 - H[np.ix_(keep, rm)] @ np.linalg.solve(H11, b1)
            U, singular, _ = np.linalg.svd((Hm + Hm.T) / 2)
            root = np.diag(np.sqrt(np.maximum(singular, 0))) @ U.T
            prior = MarginFactor([active[i] for i in keep_columns], root, np.linalg.lstsq(root.T, bm, rcond=None)[0])
            self.add_factor(prior)
        for state in self.states:
            if state.gid in remove: state.status = "Margin"
        for factor in self.factors:
            if any(state.gid in remove for state in factor.states): factor.status = "Margin"
        for index, state in enumerate(self.active_states, 1): state.lid = index
        self.win_size = len(self.active_states); return self

    def mar_measurements(self, gid):
        state = next(item for item in self.states if item.gid == gid)
        self.normal_equation()
        active = self.active_states
        size = len(state.value)
        current_index = next(index for index, item in enumerate(active) if item.gid == gid)
        columns = slice(current_index * size, (current_index + 1) * size)
        information = self.J[:, columns].T @ self.J[:, columns]
        related = [factor for factor in self.factors if factor.status != "Margin"
                   and any(item.gid == gid for item in factor.states)]
        for factor in related:
            factor.status = "Margin"
        self.add_factor(__import__("kfv_fgo.factors", fromlist=["PositionFactor"]).PositionFactor(
            [state], state.value.copy(), information + np.eye(size) * 1e-12))
        self.margin_meas_time = 0.0
        return self
