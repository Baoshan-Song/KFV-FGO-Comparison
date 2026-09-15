import numpy as np
from .Estimator import Estimator
from ..fgo.factor_graph import FactorGraph
from ..fgo.factor import State, PositionFactor, RangeFactor, PropagateFactor
from config.config import motion
from .stage_timing import StageTiming

class FgoEstimator(Estimator):
    def __init__(self, config, data, *, profile_stages=False):
        self.config = config
        self.data = data
        self.profile_stages = profile_stages

    def run(self):
        cfg = self.config.fgo
        n = self.data["num_steps"]
        timing = StageTiming(self.profile_stages, n)
        emitters = self.data["emitter_positions"]
        range_meas = self.data["toa_measurements"]

        is_imitate_kfv = getattr(cfg, "imitate_kfv", False) or getattr(cfg, "imitate_KFV", False)

        # 1. Prior Factor
        with timing.measure("add_state_factor_ms", 1):
            pos = self.data["true_positions"][:, 0]
            vel = self.data["true_velocities"][:, 0]
            x0 = np.r_[pos, vel] + cfg.err_x0

            graph = FactorGraph(cfg)

            first_state = State(1, 1, x0.copy())
            graph.add_state(first_state)

            p0_mat = cfg.p0 if isinstance(cfg.p0, np.ndarray) else np.diag(cfg.p0)
            first_factor = PositionFactor([first_state], x0.copy(), np.linalg.inv(p0_mat))
            graph.add_factor(first_factor)

            omega_r = np.linalg.inv(cfg.r) if isinstance(cfg.r, np.ndarray) else np.array([[1.0 / cfg.r]])

        # 2. main loop (2 -> N)
        for i in range(2, n + 1):
            # latest state
            # -------------------------------------------------------------
            # A. PropagateFactor
            # -------------------------------------------------------------
            with timing.measure("add_state_factor_ms", i):
                current_state = graph.states[-1]
                new_state_value = motion(current_state.value, cfg.dt, cfg.omega)
                new_state = State(i, graph.win_size + 1, new_state_value)

                graph.add_state(new_state)
                prop_factor = PropagateFactor([current_state, new_state], cfg)
                graph.add_factor(prop_factor)

            # -------------------------------------------------------------
            # B. imitate_KFV mode: schur complement
            # -------------------------------------------------------------
            if is_imitate_kfv:
                # i-1
                with timing.measure("marginalize_ms", i):
                    graph.marginalize(i - 1)

            # -------------------------------------------------------------
            # C. add all Range factors
            # -------------------------------------------------------------
            with timing.measure("add_state_factor_ms", i):
                for emitter_idx in range(emitters.shape[1]):
                    measurement = {
                        "range": range_meas[emitter_idx, i - 1],
                        "emitter": emitters[:, emitter_idx],
                        "loss_type": cfg.robust_kernel,
                        "loss_delta": cfg.robust_delta,
                        "autoDiff": getattr(cfg, "autoDiff", False),
                    }
                    range_factor = RangeFactor([new_state], measurement, omega_r)
                    graph.add_factor(range_factor)

            # -------------------------------------------------------------
            # D. nonlinear optimization (Gauss-Newton Optimization / Estimate)
            # -------------------------------------------------------------
            # Standard FGO: estimate
            with timing.measure("estimate_ms", i):
                graph.estimate()

            # -------------------------------------------------------------
            # E. Marginalization after estimation
            # -------------------------------------------------------------
            if is_imitate_kfv:
                # imitate_KFV mode: QR decomposition
                if cfg.window_size == 1:
                    with timing.measure("marginalize_ms", i):
                        graph.mar_measurements(i)
            else:
                # Standard SW-FGO: Schur complement 
                # Retire only on overflow. W >= N retains the complete graph.
                if graph.win_size > cfg.window_size:
                    with timing.measure("marginalize_ms", i):
                        graph.marginalize(i - cfg.window_size)

        # 3. output
        est_positions = np.column_stack([s.value for s in graph.states])

        result = {
            "X": est_positions,
            "debug_info": getattr(graph, "residual_norm_all", None)
        }
        if self.profile_stages:
            result["stage_timings_ms"] = timing.totals()
            result["stage_timings_by_epoch_ms"] = timing.rows
        return result



__all__ = ["FgoEstimator"]
