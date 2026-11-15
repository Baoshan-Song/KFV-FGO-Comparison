# KFV–FGO Comparison

**Open Colab, choose a comparison, and run. No installation on your computer is needed.**

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Baoshan-Song/KFV-FGO-Comparison/blob/main/KFV-FGO-Colab.ipynb)

## Three modules, two data sources

1. Run **Setup** once. Code and dependencies load automatically, with no source selection needed.
2. Run the cell for **KFV vs KFV**, **KFV vs FGO**, or **FGO vs FGO** to open its workspace.
3. Choose **Simulation** or **Real GNSS/INS**. Configure the two method cards, then click **Run comparison**. **Use saved results** is on by default.
4. Inspect the plots and export with **Save SVG**.

Both data sources are available inside every module. There is no separate real-data page or article-specific module. The Run button starts a background calculation. Stop ends after the current epoch and keeps the completed plots. The notebook stays responsive while the worker computes. Save outputs before the cloud session expires.

### Updating an already open Colab session

After opening an updated notebook, run **Setup** again, then run the comparison cell again. An existing panel keeps its old Python callbacks until it is recreated. Setup prints the loaded `main` commit and the UI revision; the current panel shows **UI 2026-11-15.4 · Colab connected** before enabling Run.

If an old panel is frozen and even Stop does not respond, use **Runtime → Restart session**, then run Setup and the comparison cell. A page refresh alone does not reload the Python code. The new panel reports startup failures directly and acknowledges Stop immediately. Its live connection also reports a disconnected or unresponsive runtime instead of silently waiting.

### Open real-data results immediately

The three bundled recordings include precomputed results for all three modules' default settings, the matched Cauchy/Tukey comparisons, and Huber SW-FGO with windows of 5 and 20 states. These are complete independent solver runs, with their settings and provenance in [the preset manifest](kfv_fgo/precomputed/manifest.json).

With **Use saved results** enabled, an exact match opens directly without parsing RINEX, aligning IMU, or solving again. The timeline, position errors, residuals and iteration diagnostics remain available. The status explicitly says **Loaded saved results**. Timing values in those results describe the original solve, not the cache load or your current computer.

Change a parameter and the panel computes a fresh result with live progress, then caches that completed run. Uncheck **Use saved results** to force a new calculation. Cache matching includes input-file contents, solver/evaluation source and the numerical configuration. It never substitutes a nearby setting or a different recording.

Fresh real-data loading reports file reading, satellite-state computation, IMU alignment and initialization separately. Parsed observations are reused in the current session. IMU intervals use binary time indexing; this preserves the original interpolated values while avoiding repeated scans of the entire recording. The worker sends numerical results through a queue. In Colab, short browser-to-kernel requests update widgets, so progress does not depend on a background asynchronous task continuing after a cell returns. Errors are displayed in the status rather than leaving a permanent loading message.

To regenerate the bundled results locally, run `python scripts/precompute_colab.py --refresh`. The script uses the same public estimator paths as a fresh UI run.

| Module | Comparisons |
|---|---|
| KFV vs KFV | EKF, iEKF, rEKF and riEKF, independently configured |
| KFV vs FGO | Any filter variant against Re-FGO or SW-FGO |
| FGO vs FGO | Two independent SW-FGO configurations, initially 5 and 20 states |

Re-FGO is an independent graph solver with two QR marginalizations. It never substitutes a precomputed filter estimate as a factor. Setting an ordinary sliding window to one does not implement the same procedure.

## Method cards

Each parameter has an icon, a short explanation and a selection or slider. Inapplicable controls are greyed out and show their effective values.

| Filter | Iterations | Kernel | Robust scale |
|---|---|---|---|
| EKF | Fixed at **1** | Fixed at **none** | Disabled |
| iEKF | Adjustable | Fixed at **none** | Disabled |
| rEKF | Fixed at **1** | Adjustable | Enabled for a robust kernel |
| riEKF | Adjustable | Adjustable | Enabled for a robust kernel |

- **Iterations** is the maximum nonlinear updates per epoch. Repetition does not add measurements.
- **Kernel** selects none, Huber, Cauchy or Tukey and operates on whitened residuals.
- **Robust scale** is dimensionless, not metres. Smaller values downweight more measurements.
- **Window** sets the states retained after an SW-FGO solve. Re-FGO fixes it at one.
- **Numerical Jacobian** is simulation-only. Off uses analytic derivatives; on uses central finite differences.
- **Re-FGO automatically follows the left filter** in KFV vs FGO. Its prior, noise, iterations and robust loss use the effective KFV configuration. The linked settings are greyed out, and no matching switch is needed. Selecting SW-FGO opens independent graph controls. The bottom FGO vs FGO module contains only SW-FGO, for comparing windows and robust settings.

For an equivalence check, choose Re-FGO in KFV vs FGO. To study iteration, compare EKF with iEKF in KFV vs KFV. To study retained history, compare different SW-FGO windows in FGO vs FGO while keeping noise and loss settings the same.

## Shared data settings

Only controls for the selected source are shown. Both methods use the same observations and initial conditions. Optional settings are grouped in collapsible panels.

| Simulation panel | Controls |
|---|---|
| Geometry and nominal noise | Number of states including the prior, range-emitter radius (m), nominal range standard deviation (m) |
| Outlier simulation | Probability, mean (m) and standard deviation (m) of the second Gaussian component |
| Initialization and repeatability | Shared X/Y initial position offsets (m), random seed |

The emitter radius controls measurement geometry, not satellite altitude. The nominal range standard deviation is squared to obtain variance. Outlier probability zero generates nominal noise only. Keep data settings and seed unchanged for a controlled comparison.

For **Real GNSS/INS**, choose bundled UrbanNav **medium**, **deep** or **harsh**. The recording loads automatically.

| Real-data panel | Controls |
|---|---|
| Process uncertainty | Position, velocity, accelerometer-bias and receiver-clock standard deviations per GNSS interval; defaults 0.30 m, 0.15 m/s, 0.01 m/s² and 100 m |
| GNSS quality gates | Optional elevation cutoff (−90° disables it) and SNR cutoff (0 disables it) |
| Custom data and time | Optional directory and GPS–UTC offset; bundled 2021 data uses 18 s |

Process standard deviations are squared to form diagonal Q. They are discrete uncertainties, not continuous IMU noise densities. A custom directory needs `f9p_navi.obs`, `brdm.rnx` and `xsens_imu.csv` or the supplied lossless `xsens_imu.npz`. Optional `gt.txt` enables reference-error plots. For a Drive data directory, mount Drive separately before selecting its path.

The real model estimates ECEF position (3), ECEF velocity (3), body accelerometer bias (3) and receiver clock bias (1). It uses recorded AHRS attitude and first-valid-sample Euler propagation per GNSS interval. It is a pseudorange/IMU estimator comparison, not carrier-phase RTK or a complete strapdown INS. Observation variances depend on elevation and SNR, following the contributed model. The simulation noise controls do not affect real data.

Truth is used only for evaluation, interpolated at overlapping UTC timestamps without extrapolation or fitted trajectory alignment. Plots use causal online states, including SW-FGO. Epoch zero is the initial prior and is excluded from accuracy statistics. Missing GNSS intervals break the plotted curves.

## Live computation and time navigation

Both simulation and real-data comparisons update while solving. The worker alternates the two methods at each epoch, and the display refreshes up to about three times per second rather than sending a whole figure for every observation. Data loading is reported before the first result appears.

- **Play / Pause** is available in the live preview. On completion, or when loading saved results, the original paper-style viewer opens automatically with the full curve and local enlargement. Its replay and slider run in the browser without a Python request for each frame. Scrubbing pauses playback, and Play at the end starts again. The original viewer is displayed directly in the Colab output, with a shorter toolbar and captions. It stays visible after completion; the live preview is hidden only after the browser confirms that the replacement figure has rendered.
- **Follow latest** tracks the newest completed pair of updates.
- Drag **Time** backwards or forwards to examine completed epochs. This turns automatic following off, so incoming results do not move your selection.
- The time label gives elapsed seconds and the epoch index, including irregular real-data intervals. **Detail (s)** controls the local time span.
- Live views show trajectory, position error and post-fit / innovation residual norms. Grey bands identify recorded GNSS gaps. Raw observation residuals are not positioning errors.
- **Stop** is cooperative and takes effect between epoch updates. Completed curves stay available, but a stopped partial run is not labeled a complete experiment.
- The finished viewer includes iteration diagnostics and exports directly. Real-data scene enlargements are optional and off by default, keeping the timeline visible.

Live curves and batch results use the same solver generators. Solver timings exclude the time spent yielding to the other method and the live viewer. The total time you wait also includes visualization and data loading, so it should not be used as a controlled solver benchmark.

## Cauchy and Tukey initialization

Their fixed-scale weights can nearly vanish when the initial state error is much larger than the measurement noise. Tukey has exactly zero weight outside its cutoff. This is a property of redescending losses, not a reason to replace their formulas with Huber or invent a nonzero weight floor.

The UI remembers a separate target scale for each kernel. Its initial suggestions are 2.385 for Cauchy and 4.685 for Tukey, following the conventional tuning constants documented by [EViews](https://help.eviews.com/content/robustreg-Background.html) and [statsmodels](https://www.statsmodels.org/stable/generated/statsmodels.robust.norms.TukeyBiweight.html). Huber retains the previous UI default of 1. These are editable starting points, not calibrated sensor noise or universal optimal settings. Applying the same cutoff of 1 to all three losses does not give them comparable rejection strength.

The method cards now expose **Gradual scale** for Cauchy and Tukey, enabled by default in the UI. The broad-to-target strategy follows the [Ceres loss-scale continuation guidance](https://ceres-solver.readthedocs.io/latest/nnls_modeling.html#lossfunctionwrapper). Our explicit starting-scale choices are twice the median absolute whitened residual for Cauchy and four times the maximum for Tukey. The latter keeps all first-step weights at least `(15/16)^2`. These choices are an initialization policy, not statistical detection thresholds.

Within a multi-update epoch, scales decrease geometrically to the selected target on the last allowed update. A single-update rEKF uses the broad scale for that update and recomputes it at the next epoch. Nominal measurement covariance is never changed, and actual scales are saved in each iteration trace. The linked Re-FGO uses the same schedule and retains both QR marginalizations. SW-FGO derives its schedule from active measurement factors.

Disable Gradual scale to reproduce the original fixed-scale estimator. It remains off by default in the core configuration and article scripts. Continuation changes the intermediate robust objective, and cannot guarantee the global optimum or that one kernel will win on every recording. A very narrow target cutoff or inadequate initialization can still reject useful observations.

## Read and export results

| View or control | Meaning |
|---|---|
| Trajectory | Full trajectory and an equal-scale local enlargement |
| Position error | Full interval and selected detail; real data offers horizontal, signed Up and 3D error |
| Residual norm | Current post-fit or pre-fit / innovation Euclidean norm (m), not the full graph objective |
| Truth measurement error | Simulation-only option within the residual view |
| Iteration path | Actual per-iteration states and residual or objective; simulation adds range contours |
| Epoch / Detail / Play | Choose and replay the enlarged interval without recalculation |
| Scene 1 / 2 enlargements | Optional medium street-region details, shared by both methods |
| Save SVG | Export the displayed vector figure |

Blue and red identify the left and right methods. Black squares mark truth, crosses mark filters, and open circles mark graphs. Matched estimates intentionally overlap. A robust estimate may have a larger raw residual norm while its position error is smaller. Interpret residuals with position error and used-satellite count. Objectives from different windows or losses are not directly comparable.

Real runs save HTML, configuration/accuracy JSON, UTC-tagged state CSVs and NPZ files in separate timestamped folders under `artifacts/colab/`. Download these from Colab's Files panel. NPZ contains causal states, predictions, covariances and per-epoch timings. SW-FGO covariance extraction is not enabled; unavailable entries are NaN.

## Article reproduction

The six-method Inside GNSS figures remain article outputs, outside the three-module UI. `python scripts/insidegnss_article.py --scene medium` retains ordinary and Huber comparisons with SW-FGO-w20, and also accepts deep or harsh. Solvers and the vector renderer remain shared with the interactive comparisons. This UI revision does not alter the article's figures or measured results.

## Supporting papers and citation

The Re-FGO theory is developed in **B. Song, R. Xu, Z. Zhan, and L.-T. Hsu, “Degeneration of sliding-window factor graph optimization into iterated extended Kalman filtering,” npj Wireless Technology, vol. 2, article 58, 2026**, [doi:10.1038/s44459-026-00069-4](https://doi.org/10.1038/s44459-026-00069-4). The GNSS/INS application builds on **W. Wen, T. Pfeifer, X. Bai, and L.-T. Hsu, “Factor graph optimization for GNSS/INS integration: A comparison with the extended Kalman filter,” NAVIGATION, vol. 68, no. 2, pp. 315–331, 2021**, [doi:10.1002/navi.421](https://doi.org/10.1002/navi.421).

Please cite the theoretical paper when using the Re-FGO comparison. If reporting results from this implementation, also identify the code version or commit and the settings used. GitHub's **Cite this repository** entry is configured through [`CITATION.cff`](CITATION.cff).

```bibtex
@article{song2026refgo,
  author = {Song, Baoshan and Xu, Ruijie and Zhan, Zhi and Hsu, Li-Ta},
  title = {Degeneration of sliding-window factor graph optimization
           into iterated extended Kalman filtering},
  journal = {npj Wireless Technology},
  volume = {2},
  number = {1},
  year = {2026},
  note = {Article 58},
  doi = {10.1038/s44459-026-00069-4}
}

@misc{song2026kfvfgo,
  author = {Song, Baoshan and {contributors}},
  title = {{KFV-FGO Comparison}: Python/Colab toolbox},
  year = {2026},
  howpublished = {\url{https://github.com/Baoshan-Song/KFV-FGO-Comparison}},
  note = {Version 0.4.0}
}
```

## Code, checks and credits

[`main`](https://github.com/Baoshan-Song/KFV-FGO-Comparison/tree/main) is the maintained Python/Colab branch. [`matlab`](https://github.com/Baoshan-Song/KFV-FGO-Comparison/tree/matlab) preserves the previous MATLAB main branch.

- Paper: [published article](https://doi.org/10.1038/s44459-026-00069-4), [v1](https://arxiv.org/html/2511.00306v1), [v2](https://arxiv.org/html/2511.00306v2). See [paper alignment](docs/paper_alignment.md) for the fixed-prediction and QR conventions.
- Thanks to **ZzPolyU** for [PR #5](https://github.com/Baoshan-Song/KFV-FGO-Comparison/pull/5), contributing GNSS/INS data readers, models, real data and saved MATLAB comparisons. Integration preserves the corrected independent graph solver. [Third-party notices](kfv_fgo/THIRD_PARTY_NOTICES.md) accompany the contribution.
- Canonical solvers live in `core/`. Real observation and propagation models live in `kfv_fgo/model/`; data adapters in `kfv_fgo/data/`. `experiment_common.py`, `real_experiment.py` and `kfv_fgo/comparison.html` drive both the GUI and article outputs.
- The [real-data verification](docs/real_validation.json) checks full-sequence state/covariance agreement for all four matched KFV/Re-FGO variants. Tests also check graph independence, affine QR priors, causal residuals, data timestamps and the contributed model probes. Historical MATLAB results remain fixtures, not substitutes for current Python outputs.

<details>
<summary>Optional local development (not required to use Colab)</summary>

Use Python 3.10 or later. In a virtual environment:

```bash
python -m pip install -e '.[real,colab,plots]' pytest
python -m unittest discover -s tests -v
python -m pytest kfv_fgo/tests -q
python scripts/check_notebook.py
python scripts/validate_paper.py --output artifacts/validation.json
python scripts/insidegnss_real.py
```

The optional Schur/window benchmark remains separate and disabled by default. It uses its legacy final-history metric rather than the GUI's online metric.

</details>

GNU General Public License v3, as specified by the original project.
