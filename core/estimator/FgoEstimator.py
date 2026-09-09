import numpy as np
from .Estimator import Estimator
from .KfvEstimator import KfvEstimator
from ..fgo.factor_graph import FactorGraph
from ..fgo.factor import State, PositionFactor, RangeFactor, PropagateFactor
from config.config import motion


class FgoEstimator(Estimator):
	def run(self):
		cfg = self.config.fgo
		count = self.data["num_steps"]
		emitters = self.data["emitter_positions"]
		if cfg.imitate_kfv and cfg.window_size == 1 and hasattr(self.config, "kfv"):
			return KfvEstimator(self.config, self.data).run()
		initial = np.r_[self.data["true_positions"][:, 0], self.data["true_velocities"][:, 0]] + cfg.err_x0
		graph = FactorGraph(cfg)
		first = State(1, 1, initial.copy())
		graph.add_state(first).add_factor(PositionFactor([first], initial.copy(), np.linalg.inv(cfg.p0)))
		for index in range(2, count + 1):
			current = graph.active_states[-1]
			new = State(index, graph.win_size + 1, motion(current.value, cfg.dt, cfg.omega))
			graph.add_state(new).add_factor(PropagateFactor([current, new], cfg))
			if cfg.window_size > 1 and len(graph.active_states) > cfg.window_size:
				graph.marginalize(graph.active_states[0].gid)
			for emitter_index in range(emitters.shape[1]):
				measurement = {
					"range": self.data["toa_measurements"][emitter_index, index - 1],
					"emitter": emitters[:, emitter_index],
					"loss_type": cfg.robust_kernel,
					"loss_delta": cfg.robust_delta,
				}
				graph.add_factor(RangeFactor([new], measurement, np.array([[1.0 / cfg.r]])))
			graph.estimate()
			if cfg.imitate_kfv and cfg.window_size == 1:
				graph.mar_measurements(new.gid)
		values = np.column_stack([state.value for state in graph.states if state.status != "Margin"])
		return {"X": values, "debug_info": graph.residual_norm_all}

__all__ = ["FgoEstimator"]
