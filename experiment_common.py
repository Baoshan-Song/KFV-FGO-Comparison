from __future__ import annotations

import json
import time
from pathlib import Path
from copy import deepcopy
import numpy as np

from config.config import Config, FilterConfig, FgoConfig, range_measurement
from core.estimator import FgoEstimator, KfvEstimator
from data.circle_eval import generate_data


def make_data(seed=7, num_steps=100, anchor_radius=200.0,
              outlier_weight=0.2, outlier_mean=0.0, outlier_sigma=10.0,
              white_sigma=0.1):
  outlier_weight = float(np.clip(outlier_weight, 0.0, 1.0))
  return generate_data(
      num_steps=num_steps,
      radius=100,
      emitter_radius=anchor_radius,
      gmm_weights=(1.0 - outlier_weight, outlier_weight),
      gmm_means=(0.0, outlier_mean),
      gmm_sigmas=(white_sigma, outlier_sigma),
      seed=seed,
  )


def make_config(err_x=100.0, err_y=-100.0, err_vx=0.0, err_vy=0.0,
                p0_diag=(50.0, 50.0, 1.0, 1.0), kfv_mode="EKF",
                robust_kernel="none", robust_delta=2.0,
                max_iteration=1, window_size=2, imitate_kfv=False,
                autodiff=False):
  initial_error = np.array([err_x, err_y, err_vx, err_vy], dtype=float)
  covariance = np.diag(np.asarray(p0_diag, dtype=float))
  kfv = FilterConfig(mode=kfv_mode, err_x0=initial_error.copy(), p0=covariance.copy(),
                     robust_kernel=robust_kernel, robust_delta=robust_delta,
                     max_iteration=max_iteration, window_size=window_size)
  fgo = FgoConfig(err_x0=initial_error.copy(), p0=covariance.copy(),
                  robust_kernel=robust_kernel,
                  robust_delta=robust_delta, max_iteration=max_iteration,
                  window_size=window_size, imitate_kfv=imitate_kfv,
                  autodiff=autodiff)
  return Config(kfv=kfv, fgo=fgo)


def run_estimator(kind, config, data):
  # Match the generated experiment's time grid without changing caller config.
  config = deepcopy(config)
  for cfg in (config.kfv, config.fgo):
    cfg.dt = float(np.asarray(data.get("dt", cfg.dt)).squeeze())
    cfg.omega = float(np.asarray(data.get("omega", cfg.omega)).squeeze())
  start_time = time.perf_counter()
  if kind == "kfv":
    result = KfvEstimator(config, data).run()
  elif kind == "fgo":
    result = FgoEstimator(config, data).run()
  else:
    raise ValueError(f"Unsupported estimator kind: {kind}")
  elapsed_ms = (time.perf_counter() - start_time) * 1000.0
  return result, elapsed_ms


def measurement_residuals(result, data, mode="postfit"):
  """Unweighted range residual z-h(x), in metres, evaluated causally.

  Epoch zero is an initial prior, not a processed measurement epoch. It is
  absent in every view. Truth noise is a diagnostic, never a solver input.
  """
  if mode == "postfit":
    estimates = result["X_online"]
  elif mode == "prefit":
    estimates = result["X_pred"]
  elif mode == "truth":
    estimates = np.vstack((data["true_positions"], data["true_velocities"]))
  else:
    raise ValueError(f"Unknown residual mode: {mode}")
  measurements, emitters = data["toa_measurements"], data["emitter_positions"]
  count = min(estimates.shape[1], measurements.shape[1])
  residuals = np.full((emitters.shape[1], count), np.nan)
  for index in range(1, count):
    for emitter_index, emitter in enumerate(emitters.T):
      residuals[emitter_index, index] = measurements[emitter_index, index] - range_measurement(estimates[:, index], emitter)
  return residuals


def metrics(result, truth, elapsed_ms):
  estimates = result["X_online"]
  count = min(estimates.shape[1], truth.shape[1])
  if count < 2:
    raise ValueError("At least one processed epoch is required")
  error = np.linalg.norm(estimates[:2, 1:count]-truth[:2, 1:count], axis=0)
  return dict(mse=float(np.mean(error**2)), rmse=float(np.sqrt(np.mean(error**2))),
              mae=float(np.mean(error)), max_error=float(np.max(error)),
              cp95=float(np.percentile(error, 95)), elapsed_ms=elapsed_ms,
              per_step_ms=elapsed_ms/(count-1), processed_epochs=count-1)


def print_statistics(results):
  print("\n" + "="*70)
  print("Online position statistics (processed epochs only):")
  for name, values in results.items():
    print(f"{name}: RMSE={values['rmse']:.3f} m | MAE={values['mae']:.3f} m | "
          f"CP95={values['cp95']:.3f} m | Max={values['max_error']:.3f} m | "
          f"Time={values['elapsed_ms']:.2f} ms ({values['per_step_ms']:.2f} ms/step)")
  print("="*70 + "\n")


def run_pair(data, first, second):
  print(f"\n🚀 [Run Estimator 1]: {first['name']} ({first['kind']})")
  first_result, first_time = run_estimator(first["kind"], first["config"], data)
  
  print(f"🚀 [Run Estimator 2]: {second['name']} ({second['kind']})")
  second_result, second_time = run_estimator(second["kind"], second["config"], data)
  
  truth = data["true_positions"]
  
  return {
      "data": data,
      "estimators": [
          {**first, "result": first_result,
           "residuals": measurement_residuals(first_result, data),
           "metrics": metrics(first_result, truth, first_time)},
          {**second, "result": second_result,
           "residuals": measurement_residuals(second_result, data),
           "metrics": metrics(second_result, truth, second_time)},
      ],
  }

def comparison_payload(experiment):
  data, estimators = experiment["data"], experiment["estimators"]
  results = [item["result"] for item in estimators]
  count = min(data["true_positions"].shape[1], *(r["X_online"].shape[1] for r in results))
  def clean(value):
    if isinstance(value, np.ndarray): return clean(value.tolist())
    if isinstance(value, np.generic): return clean(value.item())
    if isinstance(value, dict): return {k: clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)): return [clean(v) for v in value]
    if isinstance(value, (float, np.floating)): return float(value) if np.isfinite(value) else None
    return value
  return clean(dict(
    names=[item["name"] for item in estimators], frames=count,
    kinds=[item["kind"] for item in estimators],
    dt=float(np.asarray(data.get("dt", 1.)).squeeze()),
    emitters=data["emitter_positions"], measurements=data["toa_measurements"][:, :count],
    truth=data["true_positions"][:, :count],
    trajectories=[r["X_online"][:2, :count] for r in results],
    metrics=[item["metrics"] for item in estimators],
    residuals={mode: [measurement_residuals(r, data, mode)[:, :count] for r in results]
               for mode in ("postfit", "prefit", "truth")},
    diagnostics=[r["epoch_diagnostics"][:count] for r in results]))


def animate_pair(experiment, interval=70):
  from IPython.display import HTML
  payload = json.dumps(comparison_payload(experiment), allow_nan=False).replace("<", "\\u003c")
  template = (Path(__file__).parent/"kfv_fgo"/"comparison.html").read_text()
  return HTML(template.replace("__ID__", f"comparison_{id(experiment)}")
              .replace("__PAYLOAD__", payload).replace("__INTERVAL__", str(int(interval))))


def run_and_display(data, first, second):
  experiment = run_pair(data, first, second)
  print_statistics({item["name"]: item["metrics"]
                    for item in experiment["estimators"]})
  return experiment


