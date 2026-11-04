"""Run the bundled sliding-window configuration on a MAT or simulated dataset."""
from pathlib import Path
from config.init_settings_swfgo import load_settings
from data.circle_eval import load_data, generate_data
from experiment_common import run_estimator, metrics, print_statistics


def main():
    root = Path(__file__).parent
    config = load_settings()
    path = root / "data" / config.data_path
    data = load_data(path) if path.exists() else generate_data()
    result, elapsed = run_estimator("fgo", config, data)
    print_statistics({"SW-FGO": metrics(result, data["true_positions"], elapsed)})
    return result


if __name__ == "__main__":
    main()
