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

[Specify your license here, e.g., MIT, if any]
