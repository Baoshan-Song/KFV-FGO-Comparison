# KFV vs FGO Comparison

This repository provides a software framework for comparing **Kalman Filter Variants (KFV)** and **Factor Graph Optimization (FGO)** methods.

## Supported Methods

- **KFV**: Includes EKF, IEKF, REKF, RIEKF
- **Recursive FGO (Re-FGO)**: Equivalent to KFV formulations
- **Sliding Window FGO (SW-FGO)**: Fixed-lag smoother version of FGO

## Examples

- Comparison between KFV and Re-FGO:  
  `matlab/example_kfv_fgo_comparison.m`

- Recursive and Sliding Window FGO:  
  `matlab/example_sw_fgo.m`

- Simulated data generation:  
  `matlab/data/circle_eval.m`

## Usage

All scripts are in MATLAB. Run the example scripts to reproduce the benchmark results.

## License

The GICI-LIB software package is distributed under GPL v3 license. Users are freedom to modify and distribute the software as they see fit, provided that they adhere to the terms and conditions set forth in the license. This includes the ability to incorporate or use GICI-LIB with other software, whether for non-commercial or commercial purposes. However, any modifications or derivative works must also be distributed under the GPL v3 license, ensuring that the software remains free and accessible to all users.
