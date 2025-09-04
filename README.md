# README — KFV ⇄ FGO Demo GUI (App Designer .mlapp)

This App is a research/demo tool to accompany our paper “FGO MythBusters: Explaining how Kalman Filter variants achieve the same performance as FGO in navigation applications”. It helps readers reproduce comparisons between Kalman Filter Variants (KFV: EKF/IEKF/REKF/RIEKF) and Factor Graph Optimization (FGO), including Sliding-Window FGO (SW-FGO), under consistent assumptions.

The App provides:
- One-click runs of KFV vs. FGO on synthetic datasets.
- Visualizations: 2D trajectories, error time series, error CDF, and timing.
- Parameter panels to align KFV and FGO configurations.
- A dashboard with key metrics (RMSE, timing, settings).

Note: This is an academic, transparent prototype. It prioritizes clarity and reproducibility over packaging or UI polish.

---

## 1. Requirements

- MATLAB R2023b or later (R2024b recommended).
- App Designer (included with MATLAB).
- Core MATLAB toolboxes:
  - Statistics and Machine Learning Toolbox (for ECDF, etc.).
  - Optimization-related functions are implemented in code; no special optimization toolbox is required.
- OS: Windows/macOS/Linux supported by MATLAB.

Optional:
- Deep Learning Toolbox if you plan to explore auto-differentiation variants in your own extensions (not required for default runs).

---

## 2. Repository Structure (relevant to the App)

- gui/
  - FGO_KF_Simulator.mlapp — the main App.
- config/
  - init_settings_kfv_fgo_comparison.m — config script for KFV vs FGO experiments.
  - init_settings_swfgo.m — config script for SW-FGO experiments.
- core/
  - estimator/ — KfvEstimator, FgoEstimator and related estimators.
  - fgo/ — factor graph, state, and factor definitions.
  - filter/ — baseline KF utilities (if applicable).
- data/
  - circle_cv_gmm_L1.mat (example dataset) and others.

You do not need to modify these files to run the App.

---

## 3. Installation

1) Download or clone the project to a local folder (e.g., C:\... \code\).

2) Open MATLAB and add paths:
- Open the App (gui/FGO_KF_Simulator.mlapp) in App Designer, or simply run the App from MATLAB.
- The App’s startup function automatically adds the necessary project subfolders to the path. If you run code manually, ensure config/, core/, core/estimator/, core/filter/, core/fgo/, data/ are on MATLAB path.

3) Verify data and config:
- config/init_settings_kfv_fgo_comparison.m and config/init_settings_swfgo.m must exist.
- data/ must contain the referenced .mat dataset (default: circle_cv_gmm_L1.mat).

---

## 4. Quick Start

- Start MATLAB.
- Open and run the App: gui/FGO_KF_Simulator.mlapp.
- The App will:
  - Load default configurations (KFV + FGO; SW-FGO).
  - Attempt to auto-load the default dataset from config.
- If data is not found, click the folder button (📁) and pick a .mat dataset containing true_positions (and other fields used in the paper).

Recommended first run:
- Use the default “Comparison between KFV and FGO” tab.
- Keep default parameters (EKF, small window, mild thresholds).
- Click the green Run button (▶).

You will see:
- 2D trajectory plot (true vs. KFV vs. FGO).
- Error series (KFV vs. FGO).
- Error CDF (KFV vs. FGO).
- Processing time bars and dashboard stats.

---

## 5. UI Guide

Main areas:
- Top-left canvas (UIAxes): 2D trajectory and anchor points.
- Right panels:
  - ErrorAxes: error time series (KFV vs. FGO).
  - CPUAxes: simple timing/CPU-like breakdown visualization.
  - KFVAxes: error CDF comparison or FGO CDF in SW-FGO tab.
- Dashboard (bottom-right): algorithm summary, RMSE, and timing.

Left controls:
- ? button: usage tips.
- 📁 button: choose and load a .mat dataset.
- ▶ button: run the current tab’s pipeline.
- 🗑 button: clear data (simple cleanup).
- ⚙ button: save/load parameters to/from .mat.

Tabs:
1) Comparison between KFV and FGO
   - Select KF: KF/EKF/IEKF/REKF/MIEKF/RMIEKF (paper focuses on EKF/IEKF/REKF/RIEKF).
   - Robust kernel: Huber/None/Cauchy.
   - Robust delta, Window size, Max iter, Thres. iter, Init. Cov, Propa. Noise (Q diagonals).
   - Press ▶ to run:
     - KfvEstimator(cfg, data) → result1
     - Convert to FGO template → FgoEstimator(cfg, data) → result2
     - Show trajectories, errors, CDF, and timing.

2) SW FGO Simulation
   - Configure SW-FGO: window size, iteration, thresholds, robust kernel, P0, Q.
   - Press ▶ to run:
     - FgoEstimator(cfg, data) only (SW-FGO).
     - Plot trajectory + FGO CDF; dashboard shows SW-FGO stats.

---

## 6. Data Format

A minimal dataset (MAT-file) should include:
- true_positions: 2 x T ground truth trajectory.
- For KFV tab (if your filters rely on):
  - position_measurements: 2 x (T+1) or 2 x T observed positions (used by fallback KF demo).
  - pdr_positions: 2 x T control/propagation increments (optional in your own configs).
- For FGO demos:
  - toa_measurements: M x T (e.g., ranges), if range factors are used.
  - emitter_positions: 2 x M anchors/beacons.
- Other fields may be required depending on your config functions.

If unsure, start with data/circle_cv_gmm_L1.mat.

---

## 7. Reproducibility Tips

- Align KFV and FGO settings:
  - Use the same P0, Q, R, robust kernel type and delta, iteration caps, and threshold.
  - Use window size = 1 to emulate recursive filters with SW-FGO (as in the paper’s equivalence discussion).
- Initialization:
  - Ensure initial state (x0) and covariance (P0) match across KFV and FGO.
- Nonlinearity and non-Gaussian noise:
  - To test robustness, enable Huber/Cauchy and adjust delta in both KFV and FGO sides.
- Timing:
  - Timing shown is indicative and measured inside the App. For benchmarking, run multiple trials and average.

---

## 8. Save/Load Parameters

- Click ⚙ → Save Params to export:
  - config_kfv_fgo (KFV+FGO comparison settings).
  - config_swfgo (SW-FGO settings).
- Click ⚙ → Load Params to import and populate UI.

This makes it easy to share exact configurations used for figures.

---

## 9. Common Issues and Fixes

- “state not found” or FGO classes not found:
  - Ensure core/fgo, core/estimator are on the MATLAB path. The App adds them automatically at startup; if you run standalone scripts, addpath project subfolders first.
- Dataset missing fields:
  - Use the provided example data or adapt your data to include the required fields (true_positions, emitter_positions, measurements).
- No plots/update:
  - Ensure you clicked ▶ on the correct tab and that data is loaded.
- Inconsistent results:
  - Verify KFV and FGO share the same noise matrices, kernels, thresholds, and window size assumptions.

---

## 10. How This Connects to the Paper

- The App implements the “equivalence under conditions” narrative:
  - With Markov assumption, linearized models, Gaussian noise, and window size = 1, SW-FGO reproduces KFV behavior.
  - By relaxing assumptions (larger window, robust kernels, iterations), SW-FGO generalizes beyond KFV.
- The visualizations mirror the paper’s figures:
  - Trajectories (true vs. KFV vs. FGO), error CDF, and timing trends.
- Use this App to explore:
  - EKF vs. FG-EKF; IEKF vs. FG-IEKF; REKF vs. FG-REKF; RIEKF vs. FG-RIEKF.
  - The effect of window size and iterations in SW-FGO (accuracy vs. cost).

---

## 11. Citation

If you use this App in academic work, please cite our paper:

- Baoshan Song, Ruijie Xu, Li-Ta Hsu, “FGO MythBusters: Explaining how Kalman Filter variants achieve the same performance as FGO in navigation applications,” 2025.

---

## 12. Contact

For questions, bug reports, or suggestions:
- Open an issue in the project repository (if public), or
- Contact the authors listed in the paper.

Thank you for using our KFV ⇄ FGO demo App. We hope it helps you explore and understand the theoretical and practical connections between recursive filters and factor graph optimization.


## License

The software package is distributed under GPL v3 license. Users are freedom to modify and distribute the software as they see fit, provided that they adhere to the terms and conditions set forth in the license. This includes the ability to incorporate or use the comparison codes with other software, whether for non-commercial or commercial purposes. However, any modifications or derivative works must also be distributed under the GPL v3 license, ensuring that the software remains free and accessible to all users.
