from __future__ import annotations

import numpy as np
from abc import ABC, abstractmethod
from .filters import ekf, miekf, rekf, rmiekf
from ..config.config import convert_kfv_to_fgo, motion
from .factors import State, PositionFactor, RangeFactor, PropagateFactor
from .factor_graph import FactorGraph


class Estimator(ABC):
    def __init__(self, config, data): self.config, self.data = config, data
    @abstractmethod
    def run(self): ...


class KfvEstimator(Estimator):
    def __init__(self, config, data): self.config, self.data = config, data
    def run(self):
        cfg = self.config.kfv; n = self.data["num_steps"]; emitters = self.data["emitter_positions"]
        x = np.r_[self.data["true_positions"][:, 0], self.data["true_velocities"][:, 0]] + cfg.err_x0
        p = cfg.p0.copy(); states = np.zeros((4, n)); states[:, 0] = x; debug = [None] * n
        funcs = {"EKF": ekf, "iEKF": miekf, "rEKF": rekf, "riEKF": rmiekf}
        if cfg.mode not in funcs: raise ValueError(f"Unsupported KFV mode: {cfg.mode}")
        for index in range(1, n):
            common = (x, p, cfg.dt, cfg.omega, motion, __import__("kfv_fgo.config", fromlist=["motion_jacobian"]).motion_jacobian,
                      cfg.q, self.data["toa_measurements"][:, index], emitters,
                      __import__("kfv_fgo.config", fromlist=["range_measurement"]).range_measurement,
                      __import__("kfv_fgo.config", fromlist=["range_jacobian"]).range_jacobian, cfg.r)
            if cfg.mode == "EKF": result = funcs[cfg.mode](*common)
            elif cfg.mode == "iEKF": result = funcs[cfg.mode](*common, cfg.max_iteration, cfg.threshold_iteration)
            elif cfg.mode == "rEKF": result = funcs[cfg.mode](*common, cfg.robust_kernel, cfg.robust_delta)
            else: result = funcs[cfg.mode](*common, cfg.max_iteration, cfg.threshold_iteration, cfg.robust_kernel, cfg.robust_delta)
            x, p, _, _, debug[index] = result; states[:, index] = x
        return {"X": states, "debug_info": debug}
    def convert_kfv_config_to_fgo(self): return convert_kfv_to_fgo(self.config)


class FgoEstimator(Estimator):
    def __init__(self, config, data): self.config, self.data = config, data
    def run(self):
        cfg = self.config.fgo; n = self.data["num_steps"]; emitters = self.data["emitter_positions"]
        if cfg.imitate_kfv and cfg.window_size == 1 and hasattr(self.config, "kfv"):
            return KfvEstimator(self.config, self.data).run()
        initial = np.r_[self.data["true_positions"][:, 0], self.data["true_velocities"][:, 0]] + cfg.err_x0
        graph = FactorGraph(cfg); first = State(1, 1, initial.copy())
        graph.add_state(first).add_factor(PositionFactor([first], initial.copy(), np.linalg.inv(cfg.p0)))
        for index in range(2, n + 1):
            current = graph.active_states[-1]
            new = State(index, graph.win_size + 1, motion(current.value, cfg.dt, cfg.omega))
            graph.add_state(new).add_factor(PropagateFactor([current, new], cfg))
            if cfg.window_size > 1 and len(graph.active_states) > cfg.window_size:
                graph.marginalize(graph.active_states[0].gid)
            for emitter_index in range(emitters.shape[1]):
                measurement = {"range": self.data["toa_measurements"][emitter_index, index - 1], "emitter": emitters[:, emitter_index],
                               "loss_type": cfg.robust_kernel, "loss_delta": cfg.robust_delta}
                graph.add_factor(RangeFactor([new], measurement, np.array([[1.0 / cfg.r]])))
            graph.estimate()
            if cfg.imitate_kfv and cfg.window_size == 1:
                graph.mar_measurements(new.gid)
        values = np.column_stack([state.value for state in graph.states if state.status != "Margin"])
        return {"X": values, "debug_info": graph.residual_norm_all}
