"""Four KFV modes and their original FGO templates on simulation data.

Run after installing this workspace: python examples/example_kfv_fgo_simulation.py
"""

import argparse
from copy import deepcopy
from pathlib import Path

from kfv_fgo.config.settings import (
    convert_kfv_to_fgo,
    load_config,
    resolve_data_path,
    resource_path,
    workspace_root,
)
from kfv_fgo.core.estimator import FgoEstimator, KfvEstimator
from kfv_fgo.core.evaluation import save_results
from kfv_fgo.data.simulation import load_data


def run(
    config_path=None,
    output=None,
    *,
    modes=("EKF", "iEKF", "rEKF", "riEKF"),
    iterations=None,
    kernel=None,
    plots=True,
):
    """Load inputs, run the estimator, evaluate and save complete histories."""
    cfg = load_config(config_path or resource_path("kfv_fgo_comparison.json"))
    if cfg.data_mode != "sim":
        raise ValueError("This example requires simulation configuration")
    if iterations is not None:
        cfg.kfv.max_iteration = iterations
    if kernel is not None:
        cfg.kfv.robust_kernel = kernel

    data = load_data(resolve_data_path(cfg))
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
        else workspace_root() / "results" / "kfv_fgo_simulation"
    )
    save_results(output, results, simulation_truth=data["true_positions"], plots=plots)
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
    args = parser.parse_args(argv)
    return run(
        args.config,
        args.output,
        iterations=args.iterations,
        kernel=args.kernel,
        plots=not args.no_plots,
        modes=args.modes,
    )


if __name__ == "__main__":
    main()
