from __future__ import annotations

import json
import time
import numpy as np
from IPython.display import HTML

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
                max_iteration=1, window_size=1, imitate_kfv=False,
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
  start_time = time.perf_counter()
  if kind == "kfv":
    result = KfvEstimator(config, data).run()
  elif kind == "fgo":
    result = FgoEstimator(config, data).run()
  else:
    raise ValueError(f"Unsupported estimator kind: {kind}")
  elapsed_ms = (time.perf_counter() - start_time) * 1000.0
  return result, elapsed_ms


def measurement_residuals(result, data):
  estimates = result["X"]
  measurements = data["toa_measurements"]
  emitters = data["emitter_positions"]
  count = min(estimates.shape[1], measurements.shape[1])
  residuals = np.zeros((emitters.shape[1], count))
  for index in range(count):
    for emitter_index, emitter in enumerate(emitters.T):
      residuals[emitter_index, index] = (
          measurements[emitter_index, index]
          - range_measurement(estimates[:, index], emitter)
      )
  return residuals


def metrics(result, truth, elapsed_ms):
  est_x = result["X"]
  count = min(est_x.shape[1], truth.shape[1])
  
  # 🔍 [Diagnostic 1]: Print actual shape of current estimator output matrix X
  print(f"  [DEBUG Metrics] Estimated Trajectory Shape: {est_x.shape} | Ground Truth Shape: {truth.shape} | Aligned Frames Count: {count}")
  
  error = np.linalg.norm(est_x[:2, :count] - truth[:2, :count], axis=0)
  return {
      "mse": float(np.mean(error ** 2)),
      "rmse": float(np.sqrt(np.mean(error ** 2))),
      "mae": float(np.mean(error)),
      "max_error": float(np.max(error)),
      "cp95": float(np.percentile(error, 95)),
      "elapsed_ms": elapsed_ms,
      "per_step_ms": elapsed_ms / count if count > 0 else 0.0,
  }


def print_statistics(results):
  print("\n" + "="*70)
  print("Statistical Results and Dimensional Diagnostics List:")
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

def animate_pair(experiment, interval=70):
  data = experiment["data"]
  truth = data["true_positions"]
  names = [item["name"] for item in experiment["estimators"]]
  results = [item["result"] for item in experiment["estimators"]]
  residuals = [item["residuals"] for item in experiment["estimators"]]
  
  # Explicitly use maximum valid column count between ground truth and results (100 frames)
  frame_count = min(truth.shape[1], *(result["X"].shape[1] for result in results))
  metrics_list = [item["metrics"] for item in experiment["estimators"]]
  
  payload = {
      "truth": truth[:, :frame_count].tolist(),
      "trajectories": [result["X"][:2, :frame_count].tolist() for result in results], # Extract 2D positions only (x, y)
      "residuals": [values[:, :frame_count].tolist() for values in residuals],
      "names": names,
      "frames": frame_count,
      "metrics": metrics_list,
  }
  data_json = json.dumps(payload, separators=(",", ":"))
  animation_id = f"estimator_animation_{id(experiment)}"
  
  return HTML(f"""
<div id="{animation_id}" style="font-family: sans-serif; max-width: 1000px">
  <div style="display:flex; gap:12px; flex-wrap:wrap"></div>
  <div style="display:flex; align-items:center; gap:8px; margin:6px 0">
    <button class="play" style="padding:4px 12px; cursor:pointer">Play</button>
    <input class="frame" type="range" min="0" max="{frame_count - 1}" value="0" style="flex:1">
    <span class="label" style="font-weight:bold; min-width:80px; text-align:right">0 / {frame_count - 1}</span>
  </div>
  <div class="panels" style="display:grid; grid-template-columns:repeat(2,minmax(0,1fr)); gap:10px"></div>
</div>
<script>
(() => {{
  const root = document.getElementById({json.dumps(animation_id)});
  const data = {data_json};
  const colors = ["#1769aa", "#d95f02", "#1b9e77", "#7570b3", "#e7298a", "#66a61e"];
  const frameInput = root.querySelector('.frame');
  const label = root.querySelector('.label');
  const play = root.querySelector('.play');
  const panels = root.querySelector('.panels');
  
  const canvases = data.names.map((name, index) => {{
    const m = data.metrics[index];
    const panel = document.createElement('div');
    panel.innerHTML = `
      <div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:4px">
        <b style="font-size:14px; color:#222">${{name}}</b>
        <span style="font-size:11px; color:#555; background:#f0f0f0; padding:2px 6px; border-radius:4px">
          CP95: <b>${{m.cp95.toFixed(2)}} m</b> | Time: <b>${{m.elapsed_ms.toFixed(1)}} ms</b> (${{m.per_step_ms.toFixed(2)}} ms/step)
        </span>
      </div>
      <canvas width="480" height="520" style="width:100%; border:1px solid #ccc; background:#fafafa; border-radius:4px"></canvas>
    `;
    panels.appendChild(panel);
    return panel.querySelector('canvas').getContext('2d');
  }});

  // 1. Dynamically compute overall maximum coordinate bounds across all trajectories and ground truth to avoid canvas clipping
  const allX = [...data.truth[0], ...data.trajectories.flatMap(t => t[0])];
  const allY = [...data.truth[1], ...data.trajectories.flatMap(t => t[1])];
  const maxX = Math.max(...allX.map(Math.abs)), maxY = Math.max(...allY.map(Math.abs));
  const maxPosition = Math.max(maxX, maxY, 1) * 1.15; // Leave 15% padding margin

  const maxResidual = Math.max(...data.residuals.flat(2).map(Math.abs), 1) * 1.15;

  function draw(context, index, frame) {{
    const width = context.canvas.width, height = context.canvas.height;
    context.clearRect(0, 0, width, height);

    // ================= Trajectory Plot =================
    const topHeight = 310, left = 55, right = 20, bottom = 40;
    const mapX = value => left + (value + maxPosition) / (2 * maxPosition) * (width - left - right);
    const mapY = value => topHeight - bottom - (value + maxPosition) / (2 * maxPosition) * (topHeight - bottom - 15);

    context.strokeStyle = '#e0e0e0'; context.lineWidth = 1;
    context.fillStyle = '#666'; context.font = '10px sans-serif'; context.textAlign = 'center';

    const posTicks = 5;
    for (let i = 0; i <= posTicks; i++) {{
      const val = -maxPosition + (2 * maxPosition * i / posTicks);
      const px = mapX(val), py = mapY(val);

      context.beginPath(); context.moveTo(px, 15); context.lineTo(px, topHeight - bottom); context.stroke();
      context.fillText(Math.round(val), px, topHeight - bottom + 14);

      context.beginPath(); context.moveTo(left, py); context.lineTo(width - right, py); context.stroke();
      context.textAlign = 'right';
      context.fillText(Math.round(val), left - 6, py + 3);
      context.textAlign = 'center';
    }}

    context.strokeStyle = '#888'; context.lineWidth = 1.2;
    context.strokeRect(left, 15, width - left - right, topHeight - bottom - 15);

    function path(series, color, mapperX, mapperY) {{
      context.strokeStyle = color; context.lineWidth = 2; context.beginPath();
      for (let j = 0; j <= frame; j++) {{
        const x = mapperX(series[0][j], j);
        const y = mapperY(series[1][j], j);
        if (j === 0) context.moveTo(x, y); else context.lineTo(x, y);
      }}
      context.stroke();
    }}

    // Draw ground truth and current estimator trajectory
    path(data.truth, '#222', mapX, mapY);
    path(data.trajectories[index], colors[index], mapX, mapY);

    const current = data.trajectories[index];
    context.fillStyle = colors[index]; context.beginPath();
    context.arc(mapX(current[0][frame]), mapY(current[1][frame]), 5, 0, 2 * Math.PI); context.fill();

    context.fillStyle = '#111'; context.font = 'bold 11px sans-serif';
    context.fillText('X Position [m]', left + (width - left - right) / 2, topHeight - 8);
    context.save();
    context.translate(15, 15 + (topHeight - bottom - 15) / 2);
    context.rotate(-Math.PI / 2);
    context.fillText('Y Position [m]', 0, 0);
    context.restore();

    // ================= Residual Plot =================
    const residualTop = 485, residualBottom = 345;
    const mapRX = (_, j) => left + j / Math.max(data.frames - 1, 1) * (width - left - right);
    const mapRY = value => residualTop - (value + maxResidual) / (2 * maxResidual) * (residualTop - residualBottom);

    context.strokeStyle = '#e0e0e0'; context.lineWidth = 1;
    context.fillStyle = '#666'; context.font = '10px sans-serif';
    const frameStep = Math.ceil(data.frames / 5);
    for (let f = 0; f < data.frames; f += frameStep) {{
      const rx = mapRX(0, f);
      context.beginPath(); context.moveTo(rx, residualBottom); context.lineTo(rx, residualTop); context.stroke();
      context.fillText(f, rx, residualTop + 14);
    }}

    const resTicks = 4;
    for (let i = 0; i <= resTicks; i++) {{
      const rVal = -maxResidual + (2 * maxResidual * i / resTicks);
      const ry = mapRY(rVal);
      context.beginPath(); context.moveTo(left, ry); context.lineTo(width - right, ry); context.stroke();
      context.textAlign = 'right';
      context.fillText(rVal.toFixed(1), left - 6, ry + 3);
    }}

    context.strokeStyle = '#888'; context.setLineDash([4, 4]); context.beginPath();
    context.moveTo(left, mapRY(0)); context.lineTo(width - right, mapRY(0)); context.stroke();
    context.setLineDash([]);

    context.strokeRect(left, residualBottom, width - left - right, residualTop - residualBottom);

    data.residuals[index].forEach((series, seriesIndex) =>
      path([series, series], colors[seriesIndex % colors.length], mapRX, mapRY)
    );

    context.fillStyle = '#111'; context.font = 'bold 11px sans-serif'; context.textAlign = 'center';
    context.fillText('Time Step (k)', left + (width - left - right) / 2, residualTop + 28);
    context.save();
    context.translate(15, residualBottom + (residualTop - residualBottom) / 2);
    context.rotate(-Math.PI / 2);
    context.fillText('Range Residual [m]', 0, 0);
    context.restore();
  }}

  function render(frame) {{
    canvases.forEach((context, index) => draw(context, index, frame));
    label.textContent = `${{frame}} / ${{data.frames - 1}}`;
  }}

  let timer = null;
  frameInput.addEventListener('input', () => render(Number(frameInput.value)));
  play.addEventListener('click', () => {{
    if (timer) {{
      clearInterval(timer); timer = null; play.textContent = 'Play'; return;
    }}
    play.textContent = 'Pause';
    timer = setInterval(() => {{
      const next = Number(frameInput.value) + 1;
      if (next >= data.frames) {{
        clearInterval(timer); timer = null; play.textContent = 'Play'; return;
      }}
      frameInput.value = next;
      render(next);
    }}, {interval});
  }});

  render(0);
}})();
</script>
""")

def run_and_display(data, first, second):
  experiment = run_pair(data, first, second)
  print_statistics({item["name"]: item["metrics"]
                    for item in experiment["estimators"]})
  return experiment


