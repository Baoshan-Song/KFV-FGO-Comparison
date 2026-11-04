# KFV–FGO Comparison

**Open Colab, choose a comparison, and press ▶. No installation on your computer is needed.**

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Baoshan-Song/KFV-FGO-Comparison/blob/main/KFV-FGO-Colab.ipynb)

Compare Kalman filter variants (KFV), recursive factor graph optimization (Re-FGO), and sliding-window FGO on a range simulation or **real GNSS/INS measurements**. Change parameters in the form, then inspect trajectories, position errors, residuals and individual nonlinear updates in the interactive viewer. Figures export as editable SVGs.

## Start in three steps

1. Click **Open in Colab** above. Sign in to Google if Colab asks.
2. Run **Setup** once. It prepares the cloud session and loads the code automatically. Leave its defaults unchanged.
3. Run a comparison cell. For the real-data example, go directly to **4. Real GNSS/INS** and choose a preset. The bundled data loads automatically, with no upload or separate download step.

You can rerun just the comparison cell after changing its parameters. There is no need to rerun Setup for each experiment. The runtime's stop button interrupts a calculation. Save results before the Colab session expires.

## What can I compare?

| Comparison | What it reveals |
|---|---|
| EKF ↔ iEKF | Effect of repeated measurement linearization |
| EKF ↔ rEKF, iEKF ↔ riEKF | Effect of robust measurement weighting |
| KFV ↔ matched Re-FGO | Whether independent gain-based and graph-based updates retain the same information |
| Re-FGO ↔ SW-FGO | Effect of retaining and relinearizing past states |
| SW-FGO ↔ SW-FGO | Effect of window size, iteration budget or robust loss |

**KFV modes.** EKF uses one update and ordinary least squares. iEKF repeats the measurement update with a fixed prediction. rEKF uses one robust update. riEKF combines iteration and robust weighting.

**Re-FGO is an independent graph solver.** It eliminates the previous state before the measurement update, then QR-compresses the solved prior–measurement system. It never calls the filter or substitutes a precomputed filter result as a factor. Merely setting a sliding window to one is not the same algorithm.

### Real-data presets

| Preset | Left / blue | Right / red | Question |
|---|---|---|---|
| Matched EKF | EKF | Re-FGO configured as EKF | Do single-update solutions coincide? |
| Matched riEKF | riEKF | Re-FGO configured as riEKF | Does correspondence hold with iteration and robust loss? |
| EKF vs riEKF | EKF | riEKF | What changes when iteration and robust weighting are enabled together? |
| Re-FGO vs SW-FGO | Robust Re-FGO | Robust sliding-window FGO | What does re-estimating past states change? |
| Custom | Your choice | Your choice | Compare any two available methods |

The method selectors below the preset apply **only to Custom**. Other presets set their methods explicitly. The notebook prints the effective configurations in the saved JSON.

## Real GNSS/INS parameters

| Form parameter | Meaning and effect |
|---|---|
| `data_directory` | Blank uses the bundled UrbanNav excerpt. A custom directory needs `f9p_navi.obs`, `brdm.rnx` and `xsens_imu.csv` in the supplied format. Optional `gt.txt` enables reference errors. A Drive path needs Drive mounted in Setup. |
| `iterations` | Maximum nonlinear updates per epoch. Used by iEKF/riEKF and the corresponding graph mode. EKF/rEKF always perform one update. More iterations do not add new measurements. |
| `robust_kernel` | `none` keeps ordinary least squares. `huber`, `cauchy`, and `tukey` reduce the influence of large **whitened** residuals. Non-robust filter modes ignore this setting. |
| `robust_delta` | Dimensionless scale of the robust loss. A smaller value downweights more measurements. It is not a residual threshold in metres. |
| `window_epochs` | States retained **after** each SW-FGO solve. An overflow epoch is optimized before the oldest state is removed. Re-FGO always retains one. Larger windows require more computation. |
| `left_method`, `right_method` | Methods compared by the Custom preset. Each can be a KFV, Re-FGO or SW-FGO. |
| `recursive_filter` | In Custom, determines the update count and loss convention that graph methods use. For example, `riEKF` enables iterated robust graph updates. |
| `minimum_elevation_deg` | Optional satellite elevation cutoff. Default −90 disables the explicit gate. |
| `minimum_snr_dbhz` | Optional signal-strength cutoff (dB-Hz). Default 0 disables the explicit gate. Nonfinite or unavailable observations are always excluded. |
| `position_process_std_m` | Per-GNSS-interval position process standard deviation, default 0.30 m per ECEF axis. |
| `velocity_process_std_mps` | Per-interval velocity process standard deviation, default 0.15 m/s per ECEF axis. |
| `bias_process_std_mps2` | Per-interval body accelerometer-bias process standard deviation, default 0.01 m/s² per axis. |
| `clock_process_std_m` | Per-interval receiver-clock-bias process standard deviation, default 100 m. |
| `gps_minus_utc_s` | GPS–UTC offset for the recording. The bundled 2021 data uses 18 s. IMU and truth timestamps are UTC. |

The four process standard deviations are **squared** to form diagonal `Q`. They are discrete process uncertainties, not continuous IMU noise densities. Both methods receive the same `Q`, initial state and prior. JSON configuration additionally exposes `P0`, the initial offset `errX0`, convergence tolerance and the elevation/SNR weighting coefficients. Default prior standard deviations are 5 m for position, 5 m/s for velocity, 500 m/s² for accelerometer bias and 100 m for clock bias, following the contributed experiment. These broad bias priors are experimental settings, not sensor specifications.

The real model uses elevation and SNR to form a separate pseudorange variance for each observation, recomputed at each linearization. `Tref=45`, `a=30`, `A=30`, and `Fref=10` preserve the contributed weighting formula in [`GnssObservationModel`](kfv_fgo/model/gnss_ins.py). The simulation's scalar `white_sigma` does not control real GNSS weighting.

### What this real-data experiment estimates

The supplied excerpt contains **438 GNSS epochs** and 437 reference-overlap epochs. Its ten states are ECEF position (3), ECEF velocity (3), body accelerometer bias (3), and receiver clock bias (1). Pseudoranges are corrected using broadcast satellite orbits/clocks, group delay, ionosphere and troposphere models through pyRTKLIB.

Attitude comes from recorded IMU/AHRS quaternions. Propagation retains the contributor's first-valid-sample Euler discretization over each GNSS interval, with constant per-interval `Q`. It does **not** integrate all high-rate IMU samples, estimate attitude/gyro bias, or resolve carrier-phase ambiguities. This is a controlled pseudorange/IMU estimator comparison, not an RTK or full strapdown INS implementation. The four-satellite minimum is a check in this supplied observation model.

Truth is used only for evaluation. It is interpolated in ECEF at overlapping UTC timestamps, without extrapolation or fitted trajectory alignment. Errors and exported viewer trajectories use causal states from the current update, including for SW-FGO.

## Simulation parameters

Cells 1–3 retain the adjustable range/motion simulation.

| Parameter | Meaning |
|---|---|
| `steps` | Number of states including the initial prior. The circular motion rate changes with this number; both solvers use the generated rate. |
| `distance` | Radius of the fixed range-emitter arrangement (m), not a satellite altitude. Changing it changes measurement geometry. |
| `seed` | Random realization. Keep it fixed when comparing configurations. |
| `white_sigma` | Standard deviation of nominal range noise (m). Its square is the nominal measurement variance in both estimators. |
| `outlier_weight` | Probability that a measurement is drawn from the second Gaussian component. 0 gives nominal Gaussian noise only. |
| `outlier_mean`, `outlier_sigma` | Mean and standard deviation of the second component (m). |
| `initial_error_x`, `initial_error_y` | Offset added to the simulated initial position (m). |
| `*_iterations`, `*_kernel`, `*_delta` | Same iteration and loss meanings as above, configured separately for each side. |
| `match_filter_settings` | In cell 2, copies the filter's effective prior, noise, iteration and loss settings into Re-FGO. Independent right-side settings are ignored while enabled. |
| `*_imitate_kfv` | Enables two-stage Re-FGO and forces a one-state window. Disable for ordinary SW-FGO. |
| `*_window_size` | SW-FGO retention window in epochs. |
| `*_numerical_jacobian` | Uses central finite differences for simulation range derivatives. This is numerical differentiation, not automatic differentiation. |

## Read the plots and export results

Blue and red distinguish the methods. Black squares mark reference positions, crosses mark filters and open circles mark graphs. Matched solutions intentionally overlap; no artificial separation is added.

| Viewer tab/control | What it shows |
|---|---|
| **Trajectory** | Complete trajectory and an equal-scale local enlargement. Real data uses East/North coordinates. |
| **Position error** | Full interval and enlarged time interval. Real data offers horizontal magnitude, signed Up, and 3D magnitude. |
| **Residual norm** | Post-fit `‖z−h(x[k|k])‖₂` or pre-fit/innovation `‖z−h(x[k|k−1])‖₂`, in metres. Real pseudorange residuals include the estimated receiver clock and model corrections. |
| **Truth measurement error** | Simulation only. It cannot be inferred from real reference positions alone because the true receiver clock and propagation errors are unknown. |
| **Iteration path** | Actual states and residual norm or objective after each iteration, including iteration 0. Real data shows local ENU paths; synthetic emitter contours are shown only in simulation. |
| **Solver difference** | Real-data 3D separation between the two online solutions. This is not error against truth. |
| **Used satellites** | Real pseudoranges actually entering each update. A downweighted range still counts as used. |
| **Epoch / Detail / Play** | Choose the end epoch and length of the enlarged interval, or animate it. These controls do not change the calculation. |
| **Save SVG** | Export exactly the visible figure as editable vector graphics for an article. |
| **Save data** | Download the plotted data, diagnostics and real-data configurations as JSON. |

Epoch 0 is the initial prior, not a measurement update. A raw residual norm can grow when a robust method stops following an outlying range. Interpret it alongside position error and the used-satellite count. The full objective includes prior, motion and robust measurement costs. Objective values from different windows or losses are not directly comparable. The GUI distinguishes meeting the step tolerance from simply exhausting the iteration budget.

The real-data cell also writes `artifacts/real-colab/` with a standalone HTML viewer, configuration/accuracy JSON, UTC-tagged state CSVs, and NPZ files containing states, predictions and covariances. Download these through Colab's Files panel. `X_online` is causal. `X` is the graph's final stored history. SW-FGO covariance extraction is not enabled; its post-initial covariance entries are NaN, rather than fabricated uncertainty values.

## Inside GNSS companion and reproduction

The article's real-data figures use the **same presets, payload and SVG renderer as Colab**. To recreate them, choose the preset stated in a caption, set its epoch/detail controls, select the named tab, and click Save SVG. The runtime numbers printed by Colab measure that session and include diagnostic collection; they will vary with hardware.

For developers, `python scripts/insidegnss_real.py` runs all four real-data presets and saves HTML/JSON/CSV/NPZ outputs to `artifacts/insidegnss-real`. There is no separate plotting model or hidden reference alignment for article figures. The older synthetic figure script remains available as an optional simulation tool, but those figures are not the current article's results.

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
  note = {Version 0.3.0}
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
