from __future__ import annotations

import numpy as np
from scipy.linalg import qr, solve_triangular
from .factor.RangeFactor import RangeFactor
from config.config import range_measurement
from ..objective import loss_cost
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
        self.last_linear_system = None
        self.last_measurement_marginalization = None
        self.iteration_trace = []

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
        """Assemble the whitened square-root system J delta = r."""
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

    def _diagnostic(self, iteration):
        """True robust objective, distinct from the IRLS surrogate norm."""
        objective, ranges = 0., []
        for factor in self.factors:
            if factor.status == "Margin":
                continue
            if isinstance(factor, RangeFactor):
                raw = factor.z["range"] - range_measurement(
                    factor.states[0].value, factor.z["emitter"])
                objective += loss_cost(np.array([raw*np.sqrt(factor.omega[0, 0])]),
                                       factor.z["loss_type"], factor.z["loss_delta"])
                if factor.states[0] is self.active_states[-1]:
                    ranges.append(raw)
            else:
                residual = factor.evaluate().b
                objective += .5*float(residual @ residual)
        return dict(iteration=iteration,
                    state=self.active_states[-1].value.copy(),
                    measurement_residual_norm_m=float(np.linalg.norm(ranges)),
                    objective=float(objective))

    def estimate(self):
        """Solve each linearized problem directly, without squaring J's condition."""
        active = self.active_states
        if not active:
            return self
        if self.config.max_iteration < 1:
            raise ValueError("max_iteration must be positive")
        state_size = len(active[0].value)
        self.iteration_trace = [self._diagnostic(0)]
        self.converged = False
        for iteration in range(1, self.config.max_iteration + 1):
            self.normal_equation()
            if self.J is None or not len(self.r):
                break
            J, r = self.J.copy(), self.r.copy()
            anchor = np.concatenate([s.value.copy() for s in active])
            self.last_linear_system = dict(A=J, b=r, anchor=anchor,
                                          gids=tuple(s.gid for s in active))
            self.residual_norm_all.append(np.r_[active[0].value.copy(), np.linalg.norm(r)])
            delta, *_ = np.linalg.lstsq(J, r, rcond=None)
            for index, state in enumerate(active):
                state.value += delta[index*state_size:(index+1)*state_size]
            self.latest_information_matrix = J.T @ J
            self.iteration_trace.append(self._diagnostic(iteration))
            self.converged = np.linalg.norm(delta)/max(len(delta), 1) < self.config.threshold_iteration
            if self.converged:
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
        """Stage 2: QR-compress the last accepted linearized objective.

        Keep the single state. Retire its consumed prediction/measurement
        factors and retain R, Q.T b and the linearization anchor. No state
        estimate is substituted as a new observation. Freezing the last solve
        also preserves the one-iteration EKF convention without an extra
        relinearization at the updated state.
        """
        active = self.active_states
        if len(active) != 1 or active[0].gid != gid:
            raise ValueError("Second marginalization requires exactly one active state")
        system = self.last_linear_system
        if system is None or system["gids"] != (gid,):
            raise ValueError("Estimate this state before second marginalization")
        A, b, anchor = system["A"], system["b"], system["anchor"]
        q, root = qr(A, mode="economic")
        n = A.shape[1]
        if root.shape != (n, n) or np.linalg.matrix_rank(root) != n:
            raise ValueError("The compressed posterior is not full rank")
        rhs = q.T @ b
        remainder = b - q @ rhs
        prior = MarginFactor(active, root, rhs, anchor)
        self._retire(set(), list(self.factors))
        self.add_factor(prior)
        inv_root = solve_triangular(root, np.eye(n))
        self.last_measurement_marginalization = dict(
            R=root.copy(), d=rhs.copy(), anchor=anchor.copy(),
            covariance=inv_root @ inv_root.T,
            constant_cost=.5*float(remainder @ remainder), source_rows=A.shape[0])
        return self


__all__ = ["FactorGraph"]
