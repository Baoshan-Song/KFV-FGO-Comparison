from pathlib import Path
from config.init_settings_swfgo import load_settings
from data.circle_eval import load_data, generate_data
from core.estimator import FgoEstimator
from example_kfv_fgo_comparison import metrics, print_statistics


def main():
    root = Path(__file__).parent
    config = load_settings()
    path = root / "data" / config.data_path
    data = load_data(path) if path.exists() else generate_data()
    result = FgoEstimator(config, data).run()
    print_statistics("FGO Estimator", metrics(result, data["true_positions"]))


if __name__ == "__main__": main()
