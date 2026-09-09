from pathlib import Path
import numpy as np
from config.init_settings_kfv_fgo_comparison import load_settings
from data.circle_eval import load_data, generate_data
from core.estimator import KfvEstimator, FgoEstimator


def metrics(result, truth):
    count = min(result["X"].shape[1], truth.shape[1])
    error = np.linalg.norm(result["X"][:2, :count] - truth[:2, :count], axis=0)
    return {"mse": float(np.mean(error ** 2)), "rmse": float(np.sqrt(np.mean(error ** 2))),
            "mae": float(np.mean(error)), "max_error": float(np.max(error)),
            "absolute_error_95": float(np.percentile(error, 95))}


def print_statistics(title, values):
    print(f"{title} Statistics:")
    print(f"  MSE:                    {values['mse']:.6f}")
    print(f"  RMSE:                   {values['rmse']:.6f}")
    print(f"  MAE:                    {values['mae']:.6f}")
    print(f"  Max Error:              {values['max_error']:.6f}")
    print(f"  95% Absolute Error:     {values['absolute_error_95']:.6f}")


def main():
    root = Path(__file__).parent
    config = load_settings()
    path = root / "data" / config.data_path
    data = load_data(path) if path.exists() else generate_data()
    kfv = KfvEstimator(config, data).run()
    fgo_config = KfvEstimator(config, data).convert_kfv_config_to_fgo()
    fgo = FgoEstimator(fgo_config, data).run()
    print_statistics("KFV Estimator", metrics(kfv, data["true_positions"]))
    print_statistics("Re-FGO Estimator", metrics(fgo, data["true_positions"]))


if __name__ == "__main__": main()
