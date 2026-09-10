import numpy as np
from .Estimator import Estimator
from ..fgo.factor_graph import FactorGraph
from ..fgo.factor import State, PositionFactor, RangeFactor, PropagateFactor
from config.config import motion


class FgoEstimator(Estimator):
    def run(self):
        cfg = self.config.fgo
        count = self.data["num_steps"]
        emitters = self.data["emitter_positions"]
        range_meas = self.data["toa_measurements"]

        # Check if imitate_kfv mode is enabled
        is_imitate_kfv = bool(getattr(cfg, "imitate_kfv", False))

        # 1. Initialize initial state x0
        pos = self.data["true_positions"][:, 0]
        vel = self.data["true_velocities"][:, 0]
        x0_true = pos.copy() if len(pos) == 4 else np.r_[pos, vel]
        x0 = x0_true + cfg.err_x0

        graph = FactorGraph(cfg)

        # 2. Add 1st state and initial prior factor (P0)
        first_state = State(1, 1, x0.copy())
        p0_mat = cfg.p0 if isinstance(cfg.p0, np.ndarray) else np.diag(cfg.p0)
        omega_p0 = np.linalg.inv(p0_mat)
        first_factor = PositionFactor([first_state], x0.copy(), omega_p0)

        graph.add_state(first_state)
        graph.add_factor(first_factor)

        omega_r = np.linalg.inv(cfg.r) if isinstance(cfg.r, np.ndarray) else np.array([[1.0 / cfg.r]])

        # 3. Time step recursion loop (2 -> N)
        for i in range(2, count + 1):
            current_state = graph.states[i - 2]  # Previous state x_{k-1}

            # -------------------------------------------------------------
            # Step A: One-step state prediction and propagate factor injection
            # -------------------------------------------------------------
            new_state_value = motion(current_state.value, cfg.dt, cfg.omega)
            new_state = State(i, graph.win_size + 1, new_state_value)

            prop_factor = PropagateFactor([current_state, new_state], cfg)
            graph.add_state(new_state)
            graph.add_factor(prop_factor)

            # -------------------------------------------------------------
            # Step B: Sliding window marginalization logic (Branching)
            # -------------------------------------------------------------
            if is_imitate_kfv:
                # [Two-Stage Marginalization - Stage 1]: Marginalize previous state (i - 1) immediately after motion update
                graph.marginalize(i - 1)
            else:
                # Standard FGO sliding window marginalization logic
                if cfg.window_size > 1:
                    if len(graph.states) == 2:
                        graph.marginalize(i - 1)
                    if graph.win_size > cfg.window_size:
                        graph.marginalize(i - cfg.window_size)

            # -------------------------------------------------------------
            # Step C: Add Range measurement factors
            # -------------------------------------------------------------
            for j in range(emitters.shape[1]):
                measurement = {
                    "range": range_meas[j, i - 1],
                    "emitter": emitters[:, j],
                    "loss_type": cfg.robust_kernel,
                    "loss_delta": cfg.robust_delta,
                    "autoDiff": cfg.autodiff,
                }
                range_factor = RangeFactor([new_state], measurement, omega_r)
                graph.add_factor(range_factor)

            # -------------------------------------------------------------
            # Step D: Current window non-linear least squares estimation (Gauss-Newton)
            # -------------------------------------------------------------
            graph.estimate()

            # -------------------------------------------------------------
            # Step E: [Two-Stage Marginalization - Stage 2] (Triggered only in imitate_kfv mode)
            # Compress all residual factors of the current frame after measurement update to generate the prior factor for the next frame
            # -------------------------------------------------------------
            if is_imitate_kfv:
                graph.mar_measurements(i)

        # Extract all historical state estimates
        est_positions = np.column_stack([s.value for s in graph.states])

        return {
            "X": est_positions,
            "debug_info": graph.residual_norm_all
        }


__all__ = ["FgoEstimator"]