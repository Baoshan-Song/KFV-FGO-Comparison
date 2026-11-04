# GNSS/INS support

[Open the main Colab GUI](https://colab.research.google.com/github/Baoshan-Song/KFV-FGO-Comparison/blob/main/KFV-FGO-Colab.ipynb). Use **4. Real GNSS/INS**.

See the [parameter guide and comparison descriptions](../README.md). This directory contains the real-data readers, observation/IMU models and saved MATLAB references contributed in [ZzPolyU’s PR #5](https://github.com/Baoshan-Song/KFV-FGO-Comparison/pull/5). All estimators use the single canonical implementation in `core/`. Re-FGO performs actual QR marginalization and does not delegate to a filter.

[Third-party notices](THIRD_PARTY_NOTICES.md) are retained. Saved MATLAB outputs are historical fixtures, not outputs of the corrected Python solver.
