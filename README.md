## KFV v.s FGO Comparison Toolbox

[![YouTube Video Demonstration](https://img.youtube.com/vi/_W2NP7gwf7s/maxresdefault.jpg)](https://www.youtube.com/watch?v=_W2NP7gwf7s)

An interactive MATLAB research and demonstration toolbox accompanying our paper:
*"Degeneration of sliding-window factor graph optimization into iterated extended Kalman filtering"*
Baoshan Song, R. Xu, Zhi Zhan, et al.
[npj Wireless Technology](https://www.nature.com/articles/s44459-026-00069-4), 2, 58 (2026)

---

## 📅 Recent Updates

* **Aug 2026:** Published in [*npj Wireless Technology* (Volume 2, Article 58)](https://doi.org/10.1038/s44459-026-00069-4).
* **Oct 2025:** Initial preprint version ("FGO MythBusters") released on [arXiv](https://arxiv.org/abs/2511.00306v1).

---

## 1. Overview & Core Motivation

While Factor Graph Optimization (FGO) have become dominant paradigms in modern robotics and multi-sensor navigation, their theoretical relationship with traditional recursive estimators remains obscured by implementation differences. This repository provides an open, transparent, and rigorous MATLAB framework demonstrating how ** sliding window FGO (SW-FGO) algebraically and statistically degenerates into Kalman Filtering Variants (KFVs)** under consistent modeling assumptions, linearizations, and identical noise configurations.

Designed for researchers, educators, and algorithmic developers, this toolbox bridges the conceptual gap between batch factor graph optimization and recursive Bayesian filtering.

## 3. Requirements

* MATLAB R2023b or later (R2024b recommended).
* App Designer (included with standard MATLAB desktop installations).
* Core Toolboxes: Statistics and Machine Learning Toolbox. *(Custom optimization algorithms are implemented natively; no external optimization toolbox required).*

---

## 4. Repository Structure

* `gui/` — `FGO_KF_Simulator.mlapp` (The main App)
* `config/` — Configuration scripts (`init_settings_kfv_fgo_comparison.m`, `init_settings_swfgo.m`)
* `core/` — Core estimators, factor definitions, and filter utilities (`estimator/`, `fgo/`, `filter/`)
* `data/` — Example datasets (e.g., `circle_cv_gmm_L1.mat`)
* `results/` — Output trajectory plots and metrics visualizers

---

## 5. Installation & Path Setup

1. Clone or download the project to a local folder.
2. Open MATLAB and ensure all subfolders are added to your path:
   ```matlab
   addpath(genpath('path/to/KFV-FGO-Comparison'));
   savepath;
   ```

## 6. Tutorial 1: Non-GUI Command-Line Workflow

### Python port

The non-GUI MATLAB workflows are also available in Python. The GUI is intentionally
not ported. Install the runtime dependencies and run the examples from the repository root:

```bash
python3 -m pip install -e .
python3 example_kfv_fgo_comparison.py
python3 example_sw_fgo.py
python3 circle_eval.py
```

The Python files follow the original tree: configuration loaders are under
`config/`, estimator and filter adapters are under `core/`, and datasets remain
under `data/`. The examples read `config/kfv_fgo_comparison.json` and
`config/swfgo.json` through the corresponding `init_settings_*.py` functions.
`load_data` reads the existing `.mat` datasets, while `generate_data` creates
the circular test data.
The translated estimators preserve the MATLAB array layout: states are `4 x T`,
measurements are `M x T`, and emitter positions are `2 x M`.

### MATLAB-compatible ReFGO and riEKF

`python example_kfv_fgo_comparison.py` now compares riEKF with ReFGO using
`imitate_kfv=True`, `window_size=1`, 10 iterations, and Huber loss on both sides.
They share the same generated measurements and base measurement variance
`R=0.01`; the entry retains seed 7 and a 0.2 outlier fraction (mean 0, sigma 10).
IPython is only needed when explicitly requesting `animate_pair`.

**Mathematical caveat: this mode reproduces the repository's MATLAB rules;
it is mathematically inaccurate as an implementation of standard fixed-prior,
fixed-noise-scale Huber MAP estimation.** In particular:

- MATLAB `margin_factor.m` leaves the prior's A and b unchanged as the state
  moves. Python now does the same in ReFGO, omitting the fixed-mean prior
  correction used by standard iterated EKF/Gauss-Newton MAP updates.
- MATLAB `RangeFactor.m` multiplies the stored information by the current
  robust weight after every evaluation. Python now accumulates these weights
  in ReFGO, matching the existing riEKF's cumulative R inflation. Standard
  Huber IRLS would recompute weights against the unchanged base noise scale.

These compatibility behaviors are enabled only by the ReFGO estimator path.
It requires window size 1; larger windows raise a configuration error.
Ordinary SW-FGO (`imitate_kfv=False`) retains anchored priors, fixed base
measurement information, and the corrected local QR Schur elimination.
Low-level factors and marginalization also default to those standard rules.
The existing Kalman filter implementation and MATLAB sources are unchanged.

Both estimators retain the last iteration's pre-update linearization for
posterior information, as in MATLAB. ReFGO still solves its own factor graph;
it does not invoke the filter or copy filter estimates. Trajectory agreement
demonstrates compatibility with this recurrence, not correctness of the
standard Huber MAP formulation or exact nonlinear Bayesian equivalence.
Compatibility factor evaluation is intentionally stateful: extra evaluations
would advance its weight history. The normal ReFGO solve evaluates each
measurement factor once per iteration and retires it after the solve.

### Optional window benchmark and phase timing (disabled by default)

Normal comparison examples do not run a window sweep or render benchmark charts.
`FgoEstimator(..., profile_stages=False)` is the default: the phase timers do
not read a clock, allocate epoch records, or add timing fields to the result.
The existing example entry points remain ordinary estimator comparisons.

The previous plotting experiment is retained separately in
`schur_window_benchmark.py`. Running it without flags only prints help.
To deliberately repeat the experiment, use:

```bash
python schur_window_benchmark.py --run-benchmark --output ../outputs/window-benchmark
# Additional opt-ins: --plot for charts; --profile-stages for phase timing.
# To redraw previously saved measurements without running any estimators:
python schur_window_benchmark.py --plot-only ../outputs/window-benchmark/summary.json
```

The optional benchmark saves CP95 and complete estimator wall/CPU timings,
raw CSV samples, trajectories, parameters and source hashes. Windows are
interleaved across repeated runs; paired policy order alternates. `--cpu`
optionally pins only the benchmark process. Interrupted wall-time samples
are flagged and retained. `--plot` requires Matplotlib; `--profile-stages`
measures add state/factor, estimate, and marginalize separately. Pie charts
require both options. No generated charts or prior experiment outputs are
included in the repository.

### Corrected local Schur marginalization

Ordinary SW-FGO absorbs only active factors incident to a removed state.
Pivoted QR eliminates the removed columns; a second QR compresses the prior
onto the retained boundary states. Unaffected factors remain exactly once,
avoiding the upstream double counting. This is square-root Schur elimination
at the chosen linearization point, not exact nonlinear marginalization.
For the adjacent-state motion and single-state range graph, the prior is at
most 4 x 4 and does not require a full-window SVD.

Position priors use +I for the A delta = b, x += delta convention. States are
retired only on window overflow after solving; W=100 keeps the full graph in
the 100-epoch benchmark. Direct discard uses the same cleanup without a prior.
Archived state values preserve the final-history CP95 metric. The optional
benchmark uses outlier_weight=0 and at most 10 GN iterations; its execution,
plotting and phase timing remain disabled unless explicitly requested.

Run regression tests with `python -m unittest discover -s tests -v`.

### Google Colab

The complete Colab workflow is in [KFV-FGO-Colab.ipynb](KFV-FGO-Colab.ipynb).
Copy the whole project folder to Google Drive, open this notebook in Colab, and
set `PROJECT_DIR` to the copied folder if its name or location differs. The
notebook mounts Drive, copies the project to `/content`, installs dependencies,
runs compile and smoke checks, executes both examples, creates `run_all.sh`,
and saves a zip archive under `MyDrive/colab_exports/`. The MATLAB GUI is not
started in Colab.

For batch processing, cluster environments, or headless script debugging, you can bypass the graphical interface and run comparisons directly via MATLAB scripts.

1. **Open MATLAB** and navigate to your local root directory containing the repository.
2. **Configure Experiment Settings:** Edit `config/init_settings_kfv_fgo_comparison.m` or `config/init_settings_swfgo.m` to adjust state noise, initial covariance (\$P\_0\$), process noise (\$Q\$), and robust kernels.
3. **Run KFV vs. FGO Comparison Script:**
4. **Run Sliding-Window FGO Script:**
5. **Inspect Outputs:** The script automatically computes RMSE metrics, saves workspace variables, and outputs generated trajectory and error plots directly to your MATLAB figure windows.

---

## 7. Tutorial 2: Interactive GUI Workflow (`FGO_KF_Simulator`)

For visual parameter tuning and real-time comparative analysis, use the interactive App Designer interface.

1. **Launch the App:** Open MATLAB, locate `gui/FGO_KF_Simulator.mlapp`, and open it in ​**App Designer**​, or type `FGO_KF_Simulator` in the MATLAB Command Window.
2. **Select Dataset:** Click the folder button (​**📁**​) on the top-left control panel to load a `.mat` dataset (e.g., `data/circle_cv_gmm_L1.mat`).
3. **Configure Tab 1 (KFV vs. FGO Comparison):**
   * Select your desired filter variant from the dropdown (EKF, IEKF, REKF, RIEKF).
   * Toggle robust loss functions (Huber, Cauchy) and adjust tuning parameters such as window size and iteration caps.
   * Click the green Run button (​**▶**​) to execute the pipeline.
4. **Configure Tab 2 (SW FGO Simulation):**
   * Switch to the **SW FGO Simulation** tab to evaluate sliding-window optimization properties independently.
   * Adjust window lengths and maximum iterations to analyze the accuracy-versus-computational cost tradeoff.
5. **Export Parameters & Results:** Click the gear icon (​**⚙**​) to save your active experiment configuration to a `.mat` file for exact manuscript figure replication.

---

## 8. Data Format

A minimal dataset (MAT-file) should include:

* `true_positions`: \$2 \\times T\$ ground truth trajectory.
* For KFV tab: `position_measurements` or `pdr_positions`.
* For FGO demos: `toa_measurements` (\$M \\times T\$) and `emitter_positions` (\$2 \\times M\$).
* Use `data/circle_cv_gmm_L1.mat` as a default starting reference.

---

## 9. Citation

If you find this repository helpful in your academic research, please cite our paper:

```
@article{song2026degeneration,
  title={Degeneration of sliding-window factor graph optimization into iterated extended Kalman filtering},
  author={Song, Baoshan and Xu, R. and Zhan, Zhi and others},
  journal={npj Wireless Technology},
  volume={2},
  pages={58},
  year={2026},
  doi={10.1038/s44459-026-00069-4}
}
```

---

## 10. License

The software package is distributed under GPL v3 license. Users are freedom to modify and distribute the software as they see fit, provided that they adhere to the terms and conditions set forth in the license. This includes the ability to incorporate or use the comparison codes with other software, whether for non-commercial or commercial purposes. However, any modifications or derivative works must also be distributed under the GPL v3 license, ensuring that the software remains free and accessible to all users.
