"""Sliding-window FGO on simulation data.

Run after installing this workspace: python examples/example_sw_fgo_simulation.py
"""

import argparse
from pathlib import Path

from kfv_fgo.config.settings import (
    load_config,
    resolve_data_path,
    resource_path,
    workspace_root,
)
from kfv_fgo.core.estimator import FgoEstimator
from kfv_fgo.core.evaluation import save_results
from kfv_fgo.data.simulation import load_data


def run(
    config_path=None,
    output=None,
    *,
    window=None,
    iterations=None,
    kernel=None,
    plots=True,
):
    """Load inputs, run the estimator, evaluate and save complete histories."""
    cfg = load_config(config_path or resource_path("swfgo.json"))
    if cfg.data_mode != "sim":
        raise ValueError("This example requires simulation configuration")
    if iterations is not None:
        cfg.fgo.max_iteration = iterations
    if kernel is not None:
        cfg.fgo.robust_kernel = kernel
    cfg.fgo.imitate_kfv = False
    if window is not None:
        cfg.fgo.window_size = window

    data = load_data(resolve_data_path(cfg))
    results = {f"SWFGO_w{cfg.fgo.window_size}": FgoEstimator(cfg, data).run()}

    output = (
        Path(output)
        if output is not None
        else workspace_root() / "results" / "sw_fgo_simulation"
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
    parser.add_argument("--window", type=int)
    args = parser.parse_args(argv)
    return run(
        args.config,
        args.output,
        iterations=args.iterations,
        kernel=args.kernel,
        plots=not args.no_plots,
        window=args.window,
    )


if __name__ == "__main__":
    main()
