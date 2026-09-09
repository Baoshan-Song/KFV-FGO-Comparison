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
		active_factors = [factor.evaluate() for factor in self.factors if factor.status != "Margin"]
		active_states = self.active_states
		state_size = len(active_states[0].value)
		state_offsets = {state.gid: index for index, state in enumerate(active_states)}
		row_count = sum(len(factor.b) for factor in active_factors)
		J = np.zeros((row_count, len(active_states) * state_size))
		residual = np.zeros(row_count)
		row = 0
		for factor in active_factors:
			count = len(factor.b)
			residual[row:row + count] = factor.b
			for state_index, state in enumerate(factor.states):
				if state.status != "Margin":
					column = state_offsets[state.gid] * state_size
					J[row:row + count, column:column + state_size] = factor.A[
						:, state_index * state_size:(state_index + 1) * state_size]
			row += count
		self.J, self.r = J, residual
		return self

	def estimate(self):
		for _ in range(self.config.max_iteration):
			self.normal_equation()
			delta, *_ = np.linalg.lstsq(self.J, -self.r, rcond=None)
			self.latest_information_matrix = self.J.T @ self.J
			state_size = len(self.active_states[0].value)
			for index, state in enumerate(self.active_states):
				state.value += delta[index * state_size:(index + 1) * state_size]
			self.residual_norm_all.append(np.r_[self.active_states[0].value, np.linalg.norm(self.r)])
			if np.linalg.norm(delta) / len(delta) < self.config.threshold_iteration:
				break
		return self

	def marginalize(self, gids):
		remove = {gids} if np.isscalar(gids) else set(gids)
		self.normal_equation()
		active_states = self.active_states
		state_size = len(active_states[0].value)
		removed = [index for index, state in enumerate(active_states) if state.gid in remove]
		remaining = [index for index, state in enumerate(active_states) if state.gid not in remove]
		if removed and remaining:
			removed_columns = np.concatenate([np.arange(index * state_size, (index + 1) * state_size) for index in removed])
			remaining_columns = np.concatenate([np.arange(index * state_size, (index + 1) * state_size) for index in remaining])
			information = self.J.T @ self.J
			vector = self.J.T @ self.r
			H11 = information[np.ix_(removed_columns, removed_columns)]
			H12 = information[np.ix_(removed_columns, remaining_columns)]
			H21 = information[np.ix_(remaining_columns, removed_columns)]
			H22 = information[np.ix_(remaining_columns, remaining_columns)]
			marginalized_h = H22 - H21 @ np.linalg.solve(H11, H12)
			marginalized_b = vector[remaining_columns] - H21 @ np.linalg.solve(H11, vector[removed_columns])
			U, singular, _ = np.linalg.svd((marginalized_h + marginalized_h.T) / 2)
			root = np.diag(np.sqrt(np.maximum(singular, 0))) @ U.T
			prior_b = np.linalg.lstsq(root.T, marginalized_b, rcond=None)[0]
			self.add_factor(MarginFactor([active_states[index] for index in remaining], root, prior_b))
		for state in self.states:
			if state.gid in remove:
				state.status = "Margin"
		for factor in self.factors:
			if any(state.gid in remove for state in factor.states):
				factor.status = "Margin"
		for index, state in enumerate(self.active_states, 1):
			state.lid = index
		self.win_size = len(self.active_states)
		return self

	def mar_measurements(self, gid):
		state = next(item for item in self.states if item.gid == gid)
		self.normal_equation()
		active_states = self.active_states
		state_size = len(state.value)
		state_index = next(index for index, item in enumerate(active_states) if item.gid == gid)
		columns = slice(state_index * state_size, (state_index + 1) * state_size)
		information = self.J[:, columns].T @ self.J[:, columns]
		for factor in self.factors:
			if factor.status != "Margin" and any(item.gid == gid for item in factor.states):
				factor.status = "Margin"
		self.add_factor(PositionFactor([state], state.value.copy(), information + np.eye(state_size) * 1e-12))
		return self

__all__ = ["FactorGraph"]
