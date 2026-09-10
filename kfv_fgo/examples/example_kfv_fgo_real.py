"""Four KFV modes and their original FGO templates on real data.

Run after installing this workspace: python examples/example_kfv_fgo_real.py
"""

import argparse
from copy import deepcopy
from pathlib import Path

from kfv_fgo.config.settings import (
    convert_kfv_to_fgo,
    load_config,
    resource_path,
    workspace_root,
)
from kfv_fgo.core.estimator import FgoEstimator, KfvEstimator
from kfv_fgo.core.evaluation import load_ground_truth, save_results
from kfv_fgo.data.gnss_ins import GnssImuDataset


def run(
    config_path=None,
    output=None,
    *,
    modes=("EKF", "iEKF", "rEKF", "riEKF"),
    iterations=None,
    kernel=None,
    data_dir=None,
    plots=True,
):
    """Load inputs, run the estimator, evaluate and save complete histories."""
    cfg = load_config(config_path or resource_path("gnss_ins.json"))
    if cfg.data_mode != "real":
        raise ValueError("This example requires real configuration")
    if data_dir is not None:
        cfg.data_path = str(Path(data_dir).expanduser().resolve())
    if iterations is not None:
        cfg.kfv.max_iteration = iterations
    if kernel is not None:
        cfg.kfv.robust_kernel = kernel

    data = GnssImuDataset(cfg)
    truth_path = data.path / "gt.txt"
    truth = load_ground_truth(truth_path) if truth_path.is_file() else None
    results = {}
    for mode in modes:
        if mode not in ("EKF", "iEKF", "rEKF", "riEKF"):
            raise ValueError(f"Unsupported filter mode: {mode}")
        current = deepcopy(cfg)
        current.kfv.mode = mode
        results[f"KFV_{mode}"] = KfvEstimator(current, data).run()
        # The original FGO template delegates to KFV. Sliding examples solve the graph.
        results[f"ReFGO_{mode}"] = FgoEstimator(convert_kfv_to_fgo(current), data).run()

    output = (
        Path(output)
        if output is not None
        else workspace_root() / "results" / "kfv_fgo_real"
    )
    save_results(output, results, truth=truth, plots=plots)
    print(f"Saved results to {output.resolve()}")
    return results


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--iterations", type=int)
    parser.add_argument("--kernel", choices=("none", "huber", "cauchy", "tukey"))
    parser.add_argument("--no-plots", action="store_true")
    parser.add_argument(
        "--modes",
        nargs="+",
        choices=("EKF", "iEKF", "rEKF", "riEKF"),
        default=("EKF", "iEKF", "rEKF", "riEKF"),
    )
    parser.add_argument("--data-dir", type=Path)
    args = parser.parse_args(argv)
    return run(
        args.config,
        args.output,
        iterations=args.iterations,
        kernel=args.kernel,
        plots=not args.no_plots,
        modes=args.modes,
        data_dir=args.data_dir,
    )


if __name__ == "__main__":
    main()
