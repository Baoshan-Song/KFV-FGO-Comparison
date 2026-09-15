from __future__ import annotations

import numpy as np
from scipy.linalg import qr
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
        self.last_marginalization = None

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
        """Eliminate only incident factors, retaining one boundary-only prior.

        At this linearization, minimize ||Jm dm + Jb db - r|| over dm.
        A rank-revealing QR projects Jb and r onto the left nullspace of Jm.
        This is the square-root form of the (generalized) Schur complement;
        it avoids normal-equation cancellation and any full-window SVD.
        Factors not incident to the removed states are neither evaluated nor
        absorbed, so they remain in the objective exactly once.
        """
        requested = {gids} if np.isscalar(gids) else set(gids)
        active = self.active_states
        removed = [state for state in active if state.gid in requested]
        remove = {state.gid for state in removed}
        if not removed:
            return self
        affected = [factor for factor in self.factors if factor.status != "Margin"
                    and any(state.gid in remove for state in factor.states)]
        neighbor_gids = {state.gid for factor in affected for state in factor.states
                         if state.gid not in remove and state.status != "Margin"}
        boundary = [state for state in active if state.gid in neighbor_gids]
        removed_dim = sum(len(state.value) for state in removed)
        boundary_dim = sum(len(state.value) for state in boundary)
        self.last_marginalization = {"removed_dim": removed_dim,
                                    "boundary_dim": boundary_dim,
                                    "local_rows": 0, "prior_rows": 0,
                                    "factor_count": len(affected), "removed_rank": 0}
        prior = None
        if affected and boundary:
            evaluated = [factor.evaluate() for factor in affected]
            row_count = sum(np.asarray(factor.b).size for factor in evaluated)
            offsets = {}
            column = 0
            for state in removed + boundary:
                offsets[state.gid] = column
                column += len(state.value)
            local_j = np.zeros((row_count, column))
            residual = np.zeros(row_count)
            row = 0
            for factor in evaluated:
                b = np.asarray(factor.b).reshape(-1)
                residual[row:row + b.size] = b
                source_column = 0
                for state in factor.states:
                    size = len(state.value)
                    if state.gid not in offsets:
                        raise ValueError("Active factor references a retired state")
                    target = offsets[state.gid]
                    local_j[row:row + b.size, target:target + size] += factor.A[
                        :, source_column:source_column + size]
                    source_column += size
                row += b.size

            q, r, _ = qr(local_j[:, :removed_dim], mode="full", pivoting=True)
            diagonal = np.abs(np.diag(r))
            tolerance = (np.finfo(local_j.dtype).eps * max(row_count, removed_dim)
                         * (diagonal.max() if diagonal.size else 0.0))
            rank = int(np.count_nonzero(diagonal > tolerance))
            projected_j = q[:, rank:].T @ local_j[:, removed_dim:]
            projected_r = q[:, rank:].T @ residual
            self.last_marginalization.update(local_rows=row_count, removed_rank=rank)
            if projected_j.shape[0]:
                # Compress to at most boundary_dim rows. The discarded residual
                # component is a state-independent constant, not information.
                prior_q, prior_a = qr(projected_j, mode="economic")
                prior_b = prior_q.T @ projected_r
                x0 = np.concatenate([state.value.copy() for state in boundary])
                prior = MarginFactor(boundary, prior_a, prior_b, x0)
                self.last_marginalization["prior_rows"] = prior_a.shape[0]

        self._retire(remove, affected)
        if prior is not None:
            self.add_factor(prior)
        return self

    def _retire(self, remove, affected):
        """Shared cleanup for marginalization and direct discard."""
        for state in self.states:
            if state.gid in remove:
                state.status, state.lid = "Margin", 0
        for factor in affected:
            factor.status = "Margin"
        # Keep retired state values for final-trajectory evaluation, but release
        # consumed factors so active processing does not scan the factor archive.
        self.factors = [factor for factor in self.factors if factor.status != "Margin"]
        for lid, state in enumerate(self.active_states, 1):
            state.lid = lid
        self.win_size = len(self.active_states)

    def discard(self, gids):
        """Remove the same states and incident factors, without a new prior."""
        requested = {gids} if np.isscalar(gids) else set(gids)
        remove = {state.gid for state in self.active_states if state.gid in requested}
        affected = [factor for factor in self.factors if factor.status != "Margin"
                    and any(state.gid in remove for state in factor.states)]
        self._retire(remove, affected)
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
