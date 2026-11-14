# KFV–FGO Comparison

**Open Colab, choose a comparison, and run. No installation on your computer is needed.**

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Baoshan-Song/KFV-FGO-Comparison/blob/main/KFV-FGO-Colab.ipynb)

## Three modules, two data sources

1. Run **Setup** once to prepare the cloud session.
2. Run the cell for **KFV vs KFV**, **KFV vs FGO**, or **FGO vs FGO** to open its workspace.
3. Choose **Simulation** or **Real GNSS/INS**. Configure the two method cards, then click **Run comparison**.
4. Inspect the plots and export with **Save SVG** or **Save data**.

Both data sources are available inside every module. There is no separate real-data page or article-specific module. Stop calculations using Colab's cell stop button. Play / Pause in the result viewer controls replay only. Save outputs before the cloud session expires.

| Module | Comparisons |
|---|---|
| KFV vs KFV | EKF, iEKF, rEKF and riEKF, independently configured |
| KFV vs FGO | Any filter variant against Re-FGO or SW-FGO |
| FGO vs FGO | Re-FGO or SW-FGO on either side, with independent settings |

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
- **Match the filter settings**, in KFV vs FGO, copies the filter's effective prior, noise, iteration and loss settings to Re-FGO. The linked controls are disabled. Selecting SW-FGO turns the link off. Disable the link to tune Re-FGO independently.

For an equivalence check, use KFV vs FGO with the link enabled. To study iteration, compare EKF with iEKF in KFV vs KFV. To study retained history, compare Re-FGO with SW-FGO in FGO vs FGO while keeping their noise and loss settings matched.

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

Process standard deviations are squared to form diagonal Q. They are discrete uncertainties, not continuous IMU noise densities. A custom directory needs `f9p_navi.obs`, `brdm.rnx` and `xsens_imu.csv` or the supplied lossless `xsens_imu.npz`. Optional `gt.txt` enables reference-error plots. Mount Drive in Setup before selecting a Drive path.

The real model estimates ECEF position (3), ECEF velocity (3), body accelerometer bias (3) and receiver clock bias (1). It uses recorded AHRS attitude and first-valid-sample Euler propagation per GNSS interval. It is a pseudorange/IMU estimator comparison, not carrier-phase RTK or a complete strapdown INS. Observation variances depend on elevation and SNR, following the contributed model. The simulation noise controls do not affect real data.

Truth is used only for evaluation, interpolated at overlapping UTC timestamps without extrapolation or fitted trajectory alignment. Plots use causal online states, including SW-FGO. Epoch zero is the initial prior and is excluded from accuracy statistics. Missing GNSS intervals break the plotted curves.

## Read and export results

| View or control | Meaning |
|---|---|
| Trajectory | Full trajectory and an equal-scale local enlargement |
| Position error | Full interval and selected detail; real data offers horizontal, signed Up and 3D error |
| Residual norm | Current post-fit or pre-fit / innovation Euclidean norm (m), not the full graph objective |
| Truth measurement error | Simulation-only option within the residual view |
| Iteration path | Actual per-iteration states and residual or objective; simulation adds range contours |
| Solver difference | Real-data 3D separation of the two online solutions, not error against truth |
| Used satellites | Pseudoranges actually entering the real-data update; downweighted observations still count |
| Epoch / Detail / Play | Choose and replay the enlarged interval without recalculation |
| Scene 1 / 2 enlargements | Optional medium street-region details, shared by both methods |
| Save SVG / Save data | Editable vector figure / plotted data and diagnostics as JSON |

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
