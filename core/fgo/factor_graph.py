from __future__ import annotations

import numpy as np
from .factor.position_factor import PositionFactor
from .factor.margin_factor import MarginFactor


class FactorGraph:
    def __init__(self, config):
        self.config = config
        self.states = []
        self.factors = []
        self.win_size = 0
        self.J = None
        self.r = None
        self.latest_information_matrix = None
        self.ls_time = 0.0
        self.margin_time = 0.0
        self.margin_meas_time = 0.0
        self.residual_norm_all = []

    @property
    def active_states(self):
        return [state for state in self.states if state.status != "Margin"]

    def add_state(self, state):
        state.lid = len(self.active_states) + 1
        self.states.append(state)
        self.win_size += 1
        return self

    def add_factor(self, factor):
        self.factors.append(factor)
        return self

    def normal_equation(self):
        """Construct normal equations, aligned with MATLAB construct_normal_equation."""
        active_factors = [factor.evaluate() for factor in self.factors if factor.status != "Margin"]
        active_states = self.active_states
        
        if not active_states or not active_factors:
            return self

        state_size = len(active_states[0].value)
        state_offsets = {state.lid: index for index, state in enumerate(active_states)}
        
        row_count = sum(len(factor.b) for factor in active_factors)
        J = np.zeros((row_count, len(active_states) * state_size))
        residual = np.zeros(row_count)
        
        row = 0
        for factor in active_factors:
            count = len(factor.b)
            residual[row:row + count] = factor.b.flatten()
            
            for state_index, state in enumerate(factor.states):
                if state.status != "Margin":
                    column = state_offsets[state.lid] * state_size
                    J[row:row + count, column:column + state_size] = factor.A[
                        :, state_index * state_size:(state_index + 1) * state_size
                    ]
            row += count

        self.J, self.r = J, residual
        return self

    def estimate(self):
        """Nonlinear least squares optimization, aligned with MATLAB estimate_fgo."""
        active_states = self.active_states
        if not active_states:
            return self

        state_size = len(active_states[0].value)
        first_gid_state = active_states[0]

        for _ in range(self.config.max_iteration):
            self.normal_equation()
            J = self.J
            r = self.r

            if J is None or len(r) == 0:
                break

            norm_r = np.linalg.norm(r)
            self.residual_norm_all.append(np.r_[first_gid_state.value.copy(), norm_r])

            # Gauss-Newton step: H = J'*J; b = J'*r
            H = J.T @ J
            b = J.T @ r
            
            try:
                delta = np.linalg.solve(H, b)
            except np.linalg.LinAlgError:
                delta, *_ = np.linalg.lstsq(H, b, rcond=None)

            # Update active state values
            for index, state in enumerate(self.active_states):
                state.value += delta[index * state_size:(index + 1) * state_size]

            self.latest_information_matrix = H.copy()

            if np.linalg.norm(delta) / max(len(delta), 1) < self.config.threshold_iteration:
                break

        return self

    def marginalize(self, gids):
        """Marginalize specified states, aligned with MATLAB construct_marginalization."""
        remove = {gids} if np.isscalar(gids) else set(gids)
        self.normal_equation()
        
        active_states = self.active_states
        if not active_states:
            return self

        state_size = len(active_states[0].value)

        removed_indices = [index for index, state in enumerate(active_states) if state.gid in remove]
        remaining_indices = [index for index, state in enumerate(active_states) if state.gid not in remove]
        remaining_states = [state for state in active_states if state.gid not in remove]

        if removed_indices and remaining_indices:
            J1 = np.hstack([self.J[:, idx * state_size:(idx + 1) * state_size] for idx in removed_indices])
            J2 = np.hstack([self.J[:, idx * state_size:(idx + 1) * state_size] for idx in remaining_indices])
            r = self.r

            H11 = J1.T @ J1
            H12 = J1.T @ J2
            H21 = J2.T @ J1
            H22 = J2.T @ J2

            b1 = J1.T @ r
            b2 = J2.T @ r

            try:
                H11_inv_H12 = np.linalg.solve(H11, H12)
                H11_inv_b1 = np.linalg.solve(H11, b1)
            except np.linalg.LinAlgError:
                H11_inv_H12, *_ = np.linalg.lstsq(H11, H12, rcond=None)
                H11_inv_b1, *_ = np.linalg.lstsq(H11, b1, rcond=None)

            H_marg = H22 - H21 @ H11_inv_H12
            b_marg = b2 - H21 @ H11_inv_b1

            U, S_vec, _ = np.linalg.svd((H_marg + H_marg.T) / 2.0)
            S_mat = np.diag(np.sqrt(np.maximum(S_vec, 0.0)))
            J0 = S_mat @ U.T
            
            try:
                r0 = np.linalg.solve(J0.T, b_marg)
            except np.linalg.LinAlgError:
                r0, *_ = np.linalg.lstsq(J0.T, b_marg, rcond=None)

            self.add_factor(MarginFactor(remaining_states, J0, r0))

        for state in self.states:
            if state.gid in remove:
                state.status = "Margin"
                state.lid = 0

        for factor in self.factors:
            if any(state.gid in remove for state in factor.states):
                factor.status = "Margin"

        lid_counter = 1
        for state in self.states:
            if state.status != "Margin":
                state.lid = lid_counter
                lid_counter += 1

        self.win_size = len(self.active_states)
        return self

    def mar_measurements(self, gid):
        """Measurement marginalization for imitate_kfv, aligned with MATLAB construct_meas_prior_factor."""
        target_state = next((item for item in self.states if item.gid == gid), None)
        if target_state is None:
            return self

        state_size = len(target_state.value)
        active_states = self.active_states
        
        # Find target state index in active states to extract corresponding block
        target_index = next((idx for idx, item in enumerate(active_states) if item.gid == gid), 0)

        # Mark all active factors associated with this GID as 'Margin'
        for factor in self.factors:
            if factor.status != "Margin":
                if any(s.gid == gid for s in factor.states):
                    factor.status = "Margin"

        if self.latest_information_matrix is not None:
            full_info = self.latest_information_matrix.copy()
            # Safely slice the block corresponding ONLY to target_state
            start_col = target_index * state_size
            end_col = (target_index + 1) * state_size
            
            if full_info.shape[0] >= end_col:
                information_matrix = full_info[start_col:end_col, start_col:end_col]
            else:
                information_matrix = np.eye(state_size)
        else:
            information_matrix = np.eye(state_size)

        # Add single PositionFactor as prior with matched dimensions
        prior_factor = PositionFactor([target_state], target_state.value.copy(), information_matrix)
        self.add_factor(prior_factor)

        return self


__all__ = ["FactorGraph"]