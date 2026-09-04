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
