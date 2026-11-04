# KFV–FGO Comparison

Python and Google Colab toolbox for comparing Kalman filter variants, recursive factor graph optimization (Re-FGO), and sliding-window FGO.

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Baoshan-Song/KFV-FGO-Comparison/blob/main/KFV-FGO-Colab.ipynb)

- **`main`** contains the maintained Python/Colab implementation, continued from `python_colab`.
- **[`matlab`](https://github.com/Baoshan-Song/KFV-FGO-Comparison/tree/matlab)** preserves the former `main`, including its MATLAB App Designer interface.
- Paper: [v1](https://arxiv.org/html/2511.00306v1), [v2](https://arxiv.org/html/2511.00306v2), [published article](https://doi.org/10.1038/s44459-026-00069-4).

The Python implementation follows the **fixed-prediction update and two-stage marginalization in v2**. The audit and finite-iteration convention are explained in [paper alignment](docs/paper_alignment.md). It intentionally corrects the previous MATLAB-compatible recurrence rather than reproducing its numerical errors. This repository demonstrates the range/motion simulation; it does not reproduce the paper's real-world GNSS/INS experiment.

## Run in Colab

Open the badge, run **Setup**, then choose a comparison cell:

1. KFV versus KFV: EKF, iEKF, rEKF, or riEKF.
2. KFV versus FGO: **Match filter settings** is enabled by default. Re-FGO receives the same prior, noise, iteration budget and loss. Disable it to explore independent configurations.
3. FGO versus FGO: compare window size, iteration count and robust loss.

Setup clones `main` into a fresh temporary directory and installs the Python package. It never deletes an existing project. For a local or Google Drive checkout, disable `USE_GITHUB` and set `LOCAL_PATH`; enable Drive mounting only if needed.

The forms expose Gaussian/GMM noise and anchor geometry. The 105 m and 1000 m anchor radii reproduce the two geometry settings; the GMM preset uses 80% σ = 0.1 m and 20% σ = 10 m, both zero mean. Changing the epoch count also changes the circular motion rate, so the runner uses the dataset's `dt` and `omega` for both estimators. The notebook uses `white_sigma²` as nominal measurement variance. No claim is made that a new random realization reproduces a published numeric table.

## What the plots mean

All animated trajectories and reported position metrics use **online** estimates `x[k|k]`, recorded before later epochs can revise the state. Epoch 0 supplies an initial prior and is excluded from measurement residuals and accuracy metrics.

| View | Definition | Units |
|---|---|---|
| Post-fit norm | `norm(z[k] − h(x[k|k]), 2)` | m |
| Innovation norm | `norm(z[k] − h(x[k|k−1]), 2)` | m |
| Measurement-error norm | `norm(z[k] − h(x_truth[k]), 2)`, simulation diagnostic only | m |
| Iteration objective | Prior/motion costs plus measurement loss over active factors | dimensionless |
| Iteration range norm | Euclidean norm of the current epoch's unweighted range residuals | m |

The figure style follows the Inside GNSS draft: black-square truth, filter crosses, and graph open circles. The two compared methods use blue and red. The Trajectory, Position error and Residual norm tabs pair a full view with an enlarged interval. The Iteration path tab pairs the range-residual contours and update path with the convergence curve. Change the epoch and detail length to inspect another region. Save SVG exports the current figure as editable vector graphics. The iteration display includes the initial iterate and the state after every update. It reports whether the step tolerance was reached; exhausting a budget is **not** reported as convergence. Objective values from different windows/losses are not directly comparable. Undamped Gauss–Newton does not guarantee monotonic descent, and robust optimization need not reduce the unweighted range norm.

For sliding windows, `result["X"]` is the final stored trajectory (each retired state at its last update), while `result["X_online"]` is causal. `X_pred` and `epoch_diagnostics` support the UI. `covariances` contains the Gaussian or local IRLS covariance for KFV/Re-FGO; SW-FGO entries after initialization are NaN because full-window covariance extraction is not enabled.

## Local installation

Python 3.10 or later, NumPy and SciPy are required. IPython/Jupyter are optional for the animated notebook.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[colab]'
python example_kfv_fgo_comparison.py
python example_kfv_kfv_comparison.py
python example_fgo_fgo_comparison.py
jupyter notebook KFV-FGO-Colab.ipynb
```

A minimal independent equivalence check:

```python
import numpy as np
from kfv_fgo import Config, FilterConfig, convert_kfv_to_fgo
from kfv_fgo import generate_data, KfvEstimator, FgoEstimator

config = convert_kfv_to_fgo(Config(kfv=FilterConfig(mode="iEKF", max_iteration=10)))
data = generate_data(num_steps=100, emitter_radius=105, seed=7)
filtered = KfvEstimator(config, data).run()
graphed = FgoEstimator(config, data).run()
np.testing.assert_allclose(filtered["X_online"], graphed["X_online"], atol=1e-8, rtol=0)
np.testing.assert_allclose(filtered["covariances"], graphed["covariances"], atol=1e-10, rtol=1e-8)
```

The graph estimator never calls the filter. `kfv_fgo` exports the same canonical implementations used by the notebook, eliminating the obsolete duplicate implementation.

## Implementation

| Location | Role |
|---|---|
| `core/filter/` | Independent Kalman-gain updates; fixed prediction within each epoch |
| `core/fgo/factor_graph.py` | Nonlinear least squares, local state elimination, second-stage QR |
| `core/fgo/factor/` | Anchored priors, motion and range factors |
| `core/objective.py` | Fixed-scale robust loss and diagnostics |
| `experiment_common.py`, `kfv_fgo/comparison.html` | Causal evaluation and interactive plots |
| `tests/` | Algebraic, estimator and UI regression checks |
| `scripts/validate_paper.py` | Four-scene numerical comparison and machine-readable report |

Re-FGO requires `imitate_kfv=True`, `window_size=1`. Stage 1 eliminates the previous state **before** adding the current measurement. It uses square-root QR elimination algebraically equivalent to the local Schur complement. Stage 2 performs an actual thin QR of the last solved prior–measurement system and stores its triangular factor, transformed right-hand side, and anchor. Historical nonlinear factors are then retired; the current state remains.

Ordinary SW-FGO optimizes first and retires on window overflow. Only factors incident to the removed states enter the local prior, preventing double counting. A window of one alone does not enable the two-stage pipeline.

The `autodiff` JSON/config key is retained for compatibility, but it selects **central finite differences**, not automatic differentiation. The Colab form labels it as a numerical Jacobian. Analytic derivatives are the default.

## Article simulation figures

```bash
python -m pip install -e '.[plots]'
python scripts/insidegnss_simulation.py --output artifacts/insidegnss
```

This exports PDF/SVG/PNG trajectory, position-error/residual, and first-update figures, plus raw arrays and a JSON results summary. The full interval and first 15 processed epochs are specified in advance. These are newly computed simulation results, separate from any GNSS/INS field data in the article.

## Verification

```bash
python -m unittest discover -s tests -v
python scripts/validate_paper.py --output docs/validation.json
python scripts/check_notebook.py
```

[Validation results](docs/validation.json) include all four Gaussian/GMM and anchor-geometry combinations. Checks compare states, covariances and convergence traces using separate filter and graph solvers. Further tests check QR's affine term and constant cost, fixed-prior iteration, stable robust weights, factor signs, causal residuals, local marginalization, and the public import path.

The optional window benchmark remains disabled by default:

```bash
python schur_window_benchmark.py --run-benchmark --output artifacts/window-benchmark
# Add --plot for figures and --profile-stages for phase timers.
```

That benchmark retains its legacy final-history accuracy metric; it is separate from the Colab online metrics. Timings include the diagnostic bookkeeping currently enabled by each estimator and are not a reproduction of published timing tables.

## License

GNU General Public License v3, as specified by the original project. Modified and redistributed versions remain subject to its terms.
