"""Original KFV shortcut for the template; sliding FGO solves the factor graph."""

from abc import ABC, abstractmethod
from time import perf_counter

import numpy as np

from ..config.settings import _config_to_dict, convert_kfv_to_fgo, validate_config
from ..model.simulation import SimulationModel
from .fgo import (
    FactorGraph,
    GnssPseudorangeFactor,
    PositionFactor,
    PropagateFactor,
    RangeFactor,
    State,
)
from .filter import filter_step


class Estimator(ABC):
    def __init__(self, config, data, *, trace=False):
        self.config, self.data, self.trace = config, data, trace
        self.real = hasattr(data, "get_measurement")

    def setup(self, cfg):
        if self.real:
            initial = self.data.initial_state + cfg.err_x0
            times = self.data.timestamps

            def propagate(x, i):
                return self.data.propagate(x, i, cfg)

            observe = self.data.linearize
        else:
            initial = (
                np.r_[
                    self.data["true_positions"][:, 0],
                    self.data["true_velocities"][:, 0],
                ]
                + cfg.err_x0
            )
            times = np.asarray(
                self.data.get("timestamps", np.arange(self.data["num_steps"]) * cfg.dt)
            )
            model = SimulationModel(cfg, self.data)
            propagate = model.propagate
            observe = model.linearize
        validate_config(self.config, initial, cfg)
        if len(times) < 1 or np.any(np.diff(times) <= 0):
            raise ValueError("Invalid dataset timestamps")
        return initial, times, propagate, observe

    def result(self, states, times, debug, started, **extra):
        return {
            "X": states,
            "timestamps": times.copy(),
            "debug_info": debug,
            "runtime_seconds": perf_counter() - started,
            "config": _config_to_dict(self.config),
            "solver": type(self).__name__,
            "backend": self.data.backend.metadata
            if self.real
            else {"model": "simulation"},
            **extra,
        }

    @abstractmethod
    def run(self):
        pass


class KfvEstimator(Estimator):
    def run(self):
        started = perf_counter()
        cfg = self.config.kfv
        x, times, propagate, observe = self.setup(cfg)
        p = cfg.p0.copy()
        values = np.empty((len(x), len(times)))
        values[:, 0] = x
        debug = [None] * len(times)
        for i in range(1, len(times)):
            try:
                prediction, F, Q = propagate(x, i)
                x, p, _, _, debug[i] = filter_step(
                    x,
                    p,
                    prediction,
                    F,
                    Q,
                    lambda state, epoch=i: observe(state, epoch),
                    mode=cfg.mode,
                    max_iter=cfg.max_iteration,
                    threshold=cfg.threshold_iteration,
                    loss_type=cfg.robust_kernel,
                    delta=cfg.robust_delta,
                    dynamic_observation=self.real,
                    trace=self.trace,
                )
            except (ValueError, np.linalg.LinAlgError, FloatingPointError) as exc:
                raise ValueError(f"KFV epoch {i + 1} at {times[i]:.6f}: {exc}") from exc
            values[:, i] = x
        return self.result(values, times, debug, started)

    def convert_kfv_config_to_fgo(self):
        return convert_kfv_to_fgo(self.config)

    convert_KFV_config_to_FGO = convert_kfv_config_to_fgo


class FgoEstimator(Estimator):
    def run(self):
        started = perf_counter()
        cfg = self.config.fgo
        # Restore the original Python template exactly, as requested. This path
        # deliberately delegates to KFV; SWFGO (imitate_kfv=False) still solves
        # factors, normal equations and marginalization below.
        if cfg.imitate_kfv and cfg.window_size == 1 and hasattr(self.config, "kfv"):
            return KfvEstimator(self.config, self.data).run()
        initial, times, propagate, observe = self.setup(cfg)
        graph = FactorGraph(cfg, trace=self.trace)
        first = State(1, 1, initial.copy())
        graph.add_state(first).add_factor(
            PositionFactor([first], initial.copy(), np.linalg.inv(cfg.p0))
        )
        debug = [None] * len(times)
        for index in range(1, len(times)):
            try:
                old = graph.states[-1]
                prediction, _, _ = propagate(old.value, index)
                new = State(index + 1, graph.win_size + 1, prediction)
                graph.add_state(new).add_factor(
                    PropagateFactor([old, new], cfg, lambda x, i=index: propagate(x, i))
                )
                # Preserve initial anchoring and pre-measurement marginalization.
                if cfg.window_size > 1 and len(graph.states) == 2:
                    graph.marginalize(old.gid)
                if graph.win_size > cfg.window_size:
                    graph.marginalize(new.gid - cfg.window_size)
                if self.real:
                    graph.add_factor(
                        GnssPseudorangeFactor(
                            [new], lambda x, i=index: observe(x, i), cfg
                        )
                    )
                else:
                    for j, emitter in enumerate(self.data["emitter_positions"].T):
                        measurement = {
                            "range": self.data["toa_measurements"][j, index],
                            "emitter": emitter,
                            "loss_type": cfg.robust_kernel,
                            "loss_delta": cfg.robust_delta,
                        }
                        graph.add_factor(
                            RangeFactor([new], measurement, np.array([[1 / cfg.r]]))
                        )
                graph.estimate()
                debug[index] = {
                    "iterations": graph.iterations,
                    "iteration_count": graph.iteration_count,
                    "residual_norm_all": graph.residual_norm_all,
                    "ls_time": graph.ls_time,
                    "margin_time": graph.margin_time,
                    "active_states": graph.win_size,
                }
                if cfg.imitate_kfv and cfg.window_size == 1:
                    graph.mar_measurements(new.gid)
            except (ValueError, np.linalg.LinAlgError, FloatingPointError) as exc:
                raise ValueError(
                    f"FGO epoch {index + 1} at {times[index]:.6f}: {exc}"
                ) from exc
        values = np.column_stack([s.value for s in graph.states])
        return self.result(
            values,
            times,
            debug,
            started,
            trajectory_semantics="State at marginalization; remaining window at final solve",
        )
