"""Optional Schur/discard benchmark; normal examples never invoke this module.

Use --run-benchmark to collect measurements, --plot to also render charts,
and --profile-stages to collect phase timings. All are disabled by default.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib
import json
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from time import perf_counter, process_time
from unittest.mock import patch

# Only an explicitly launched benchmark controls its process's BLAS threads.
# Importing helpers from tests or notebooks must not change their environment.
if __name__ == "__main__" and "--run-benchmark" in sys.argv:
    for variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
                     "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[variable] = "1"

import numpy as np
from config.config import Config, FgoConfig
from core.fgo.factor_graph import FactorGraph
from data.circle_eval import generate_data

ESTIMATOR = importlib.import_module("core.estimator.FgoEstimator")
POLICIES = ("schur", "discard")
LABELS = {"schur": "Schur complement", "discard": "Direct discard"}
COLORS = {"schur": "#2563A6", "discard": "#D88C38"}


class ComparisonGraph(FactorGraph):
    """Switch retention policy while sharing estimation and factor cleanup."""

    def __init__(self, config, policy):
        super().__init__(config)
        self.policy = policy
        self.removal_calls = 0
        self.max_solve_states = 0
        self.max_boundary_dim = 0
        self.max_local_rows = 0

    def estimate(self):
        self.max_solve_states = max(self.max_solve_states, self.win_size)
        return super().estimate()

    def marginalize(self, gids):
        self.removal_calls += 1
        if self.policy == "schur":
            super().marginalize(gids)
            stats = self.last_marginalization
            if stats is not None:
                self.max_boundary_dim = max(self.max_boundary_dim, stats["boundary_dim"])
                self.max_local_rows = max(self.max_local_rows, stats["local_rows"])
        else:
            super().discard(gids)
        return self


def run_case(data, window, policy, profile_stages=False):
    config = Config(fgo=FgoConfig(window_size=window, imitate_kfv=False,
                                 robust_kernel="none", max_iteration=10,
                                 autodiff=False))
    graphs = []

    def factory(cfg):
        graph = ComparisonGraph(cfg, policy)
        graphs.append(graph)
        return graph

    with patch.object(ESTIMATOR, "FactorGraph", factory):
        started_utc = datetime.now(timezone.utc).isoformat(timespec="milliseconds")
        cpu_started = process_time()
        started = perf_counter()
        result = ESTIMATOR.FgoEstimator(config, data, profile_stages=profile_stages).run()
        elapsed_ms = (perf_counter() - started) * 1000
        cpu_ms = (process_time() - cpu_started) * 1000
    trajectory = result["X"]
    if trajectory.shape != (4, data["num_steps"]) or not np.isfinite(trajectory).all():
        raise ValueError(f"Invalid complete trajectory: W={window}, {policy}")
    error = np.linalg.norm(trajectory[:2] - data["true_positions"], axis=0)
    graph = graphs[0]
    record = {"window": window, "policy": policy,
            "cp95_m": float(np.percentile(error, 95)), "elapsed_ms": elapsed_ms,
            "cpu_ms": cpu_ms, "started_utc": started_utc,
            "wall_cpu_ratio": elapsed_ms / cpu_ms if cpu_ms > 0 else None,
            "timing_suspect": elapsed_ms - cpu_ms > 1000 and elapsed_ms > 2 * cpu_ms,
            "removal_calls": graph.removal_calls,
            "max_solve_states": graph.max_solve_states,
            "final_active_states": graph.win_size,
            "total_gn_iterations": len(graph.residual_norm_all),
            "max_boundary_dim": graph.max_boundary_dim,
            "max_local_rows": graph.max_local_rows}
    if profile_stages:
        record.update(result["stage_timings_ms"])
        record["other_ms"] = elapsed_ms - sum(result["stage_timings_ms"].values())
        if record["other_ms"] < 0:
            raise RuntimeError("Phase timers overlap or exceed the complete estimator time")
        record["_stage_rows"] = result["stage_timings_by_epoch_ms"]
    return record, trajectory


def write_csv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def source_hashes(root):
    paths = sorted([*root.joinpath("core").rglob("*.py"),
                    *root.joinpath("config").rglob("*.py"),
                    root / "data/circle_eval.py", Path(__file__).resolve(),
                    root / "stage_timing_plots.py"])
    return {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in paths}


def set_benchmark_cpu(cpu):
    """Pin only this benchmark process; leave system power/priority unchanged."""
    if cpu is None:
        return
    if os.name == "nt":
        import ctypes
        from ctypes import wintypes
        kernel = ctypes.WinDLL("kernel32", use_last_error=True)
        kernel.GetCurrentProcess.restype = wintypes.HANDLE
        kernel.SetProcessAffinityMask.argtypes = (wintypes.HANDLE, ctypes.c_size_t)
        kernel.SetProcessAffinityMask.restype = wintypes.BOOL
        if not kernel.SetProcessAffinityMask(kernel.GetCurrentProcess(), 1 << cpu):
            raise ctypes.WinError(ctypes.get_last_error())
    elif hasattr(os, "sched_setaffinity"):
        os.sched_setaffinity(0, {cpu})
    else:
        raise RuntimeError("CPU affinity is unsupported on this platform; omit --cpu")


def render_plots(output, summary):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                         "axes.labelsize": 12, "svg.fonttype": "none"})
    windows = summary["windows"]
    lookup = {(r["window"], r["policy"]): r for r in summary["aggregates"]}
    positions = np.arange(len(windows))
    charts = [
        ("cp95_m", "Final-trajectory accuracy", "CP95 position error (m)", "cp95_by_window"),
        ("elapsed_ms", "Computation time", "Elapsed time for 100 epochs (ms)", "runtime_by_window"),
    ]
    if all("cpu_ms" in row for row in summary["aggregates"]):
        charts.append(("cpu_ms", "CPU computation time", "Process CPU time for 100 epochs (ms)", "cpu_runtime_by_window"))
    for metric, title, ylabel, name in charts:
        fig, ax = plt.subplots(figsize=(12, 6.7))
        fig.subplots_adjust(left=.10, right=.98, bottom=.23, top=.76)
        positive = [r[metric] for r in summary["aggregates"] if r[metric] > 0]
        log_scale = max(positive) / min(positive) > 100
        if log_scale:
            ax.set_yscale("log")
        upper = []
        for offset, policy in ((-.20, "schur"), (.20, "discard")):
            rows = [lookup[(w, policy)] for w in windows]
            heights = np.array([r[metric] for r in rows])
            errorbar = None
            if metric in ("elapsed_ms", "cpu_ms"):
                prefix = metric.removesuffix("_ms")
                low = np.array([r[f"{prefix}_q25_ms"] for r in rows])
                high = np.array([r[f"{prefix}_q75_ms"] for r in rows])
                errorbar = np.vstack((heights - low, high - heights))
                upper.extend(high)
            else:
                upper.extend(heights)
            bars = ax.bar(positions + offset, heights, width=.35,
                          label=LABELS[policy], color=COLORS[policy],
                          yerr=errorbar, capsize=3, zorder=3,
                          error_kw={"elinewidth": 1, "ecolor": "#263445"})
            labels = [f"{v:.2g}" if log_scale else
                      (f"{v:.3f}" if v < 1 else f"{v:.2f}" if v < 100 else f"{v:.0f}") for v in heights]
            ax.bar_label(bars, labels=labels, padding=5, fontsize=9, color="#243244")
        ax.set_xticks(positions, [str(w) for w in windows])
        ax.set_xlabel("Configured sliding-window size W (states)", labelpad=12)
        ax.set_ylabel(ylabel + (" — log scale" if log_scale else ""), labelpad=10)
        ax.grid(axis="y", alpha=.18, zorder=0)
        ax.spines[["top", "right"]].set_visible(False)
        ax.spines[["left", "bottom"]].set_color("#BDC6D1")
        ax.set_ylim((min(positive) / 2 if log_scale else 0),
                    max(upper) * (3 if log_scale else 1.22))
        fig.text(.10, .925, title, fontsize=21, fontweight="bold", color="#16283E")
        fig.text(.10, .872, "Schur complement vs. direct discard  |  Lower is better",
                 fontsize=12, color="#526173")
        ax.legend(loc="lower left", bbox_to_anchor=(0, 1.015), ncol=2,
                  frameon=False, borderaxespad=0, fontsize=11)
        detail = ("Final-history 2D error, 95th percentile; one fixed dataset (seed 7)."
                  if metric == "cp95_m" else
                  f"Median of {summary['repeats']} runs after {summary['warmups']} warm-up(s) per case; "
                  "error bars: 25th–75th percentiles. BLAS: 1 thread.")
        if metric == "cpu_ms":
            detail = (f"Process CPU time excludes pauses and time not scheduled. "
                      f"Median and 25th–75th percentiles of {summary['repeats']} runs.")
        if metric in ("elapsed_ms", "cpu_ms") and summary.get("logical_cpu") is not None:
            detail += f" CPU: {summary['logical_cpu']}."
        fig.text(.10, .107, detail, fontsize=9.5, color="#526173")
        fig.text(.10, .068,
                 ("W=100: identical batch computation in both policies; any gap reflects timing variability."
                  if metric in ("elapsed_ms", "cpu_ms") and summary.get("marginalization_algorithm") else
                  "Local Schur prior; W=100 keeps all states (batch). Window cleanup follows each solve."
                  if summary.get("marginalization_algorithm") else
                  "Source algorithm retained: initial-state removal and post-solve cleanup; W=100 is not strict batch."),
                 fontsize=9.2, color="#526173")
        fig.text(.10, .029,
                 f"100 epochs · 6 anchors · max {summary['max_iteration']} GN iterations/update · "
                 f"outlier weight {summary['data']['outlier_weight']:g} · no robust kernel",
                 fontsize=9.2, color="#526173")
        fig.savefig(output / f"{name}.png", dpi=180, facecolor="white")
        fig.savefig(output / f"{name}.svg", facecolor="white")
        plt.close(fig)
    if summary.get("profile_stages"):
        from stage_timing_plots import render_stage_plots
        render_stage_plots(output, summary)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-benchmark", action="store_true",
                        help="Explicitly enable benchmark execution and result files")
    parser.add_argument("--plot", action="store_true", help="Also render charts after the benchmark")
    parser.add_argument("--windows", type=int, nargs="+", default=[2, 5, 10, 20, 30, 50, 75, 100])
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--cpu", type=int, help="Pin benchmark to one logical CPU to reduce scheduling variation")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--profile-stages", action="store_true",
                        help="Measure three phases; pie charts additionally require --plot")
    parser.add_argument("--plot-only", type=Path, help="Re-render an existing summary.json")
    args = parser.parse_args(argv)
    if args.plot_only:
        render_plots(args.plot_only.parent,
                     json.loads(args.plot_only.read_text(encoding="utf-8")))
        return
    if not args.run_benchmark:
        parser.print_help()
        return
    if (args.repeats < 1 or args.warmups < 0 or
            len(set(args.windows)) != len(args.windows) or
            any(w < 2 or w > 100 for w in args.windows)):
        parser.error("Use unique windows in 2..100, repeats >= 1, warmups >= 0")
    if args.cpu is not None and not 0 <= args.cpu < min(os.cpu_count() or 1, 64):
        parser.error("--cpu must identify an available logical CPU below 64")
    set_benchmark_cpu(args.cpu)
    root = Path(__file__).resolve().parent
    output = args.output or (root.parent / "outputs" /
                            datetime.now().strftime("schur-window-benchmark-%Y%m%d-%H%M%S"))
    output.mkdir(parents=True, exist_ok=False)
    before = source_hashes(root)
    # Same geometry and seed as the earlier comparison, with no outliers.
    data = generate_data(num_steps=100, radius=100, emitter_radius=105,
                         gmm_weights=(1., 0.), gmm_means=(0., 0.),
                         gmm_sigmas=(.1, 10.), seed=7)
    np.savez_compressed(output / "input_data.npz", **data)
    records, aggregates, trajectories = [], [], {}
    epoch_records = []
    cp95_reference = {}
    windows = sorted(args.windows)
    for warmup in range(args.warmups):
        print(f"Warm-up {warmup + 1}/{args.warmups}: all windows", flush=True)
        for window in (windows if warmup % 2 == 0 else windows[::-1]):
            for policy in POLICIES:
                run_case(data, window, policy, args.profile_stages)
    # Spread each configuration over the complete run instead of confounding
    # window size with desktop load, CPU frequency or thermal drift over time.
    order_rng = np.random.default_rng(730)
    for repeat in range(args.repeats):
        order = order_rng.permutation(windows).tolist()
        print(f"Repeat {repeat + 1}/{args.repeats}: windows {order}", flush=True)
        for window in order:
            for policy in (POLICIES if repeat % 2 == 0 else POLICIES[::-1]):
                row, trajectory = run_case(data, window, policy, args.profile_stages)
                for epoch_row in row.pop("_stage_rows", []):
                    epoch_records.append({"window": window, "policy": policy,
                                          "repeat": repeat + 1, **epoch_row})
                key = f"{policy}_w{window}"
                if key in cp95_reference and not np.isclose(
                        row["cp95_m"], cp95_reference[key], rtol=1e-10, atol=1e-10):
                    raise RuntimeError(f"Unexpected nondeterministic CP95: {key}")
                cp95_reference[key] = row["cp95_m"]
                trajectories[key] = trajectory
                records.append({"repeat": repeat + 1, **row})
        write_csv(output / "raw_measurements.csv", records)
        if epoch_records:
            write_csv(output / "stage_timings_by_epoch.csv", epoch_records)
    for window in windows:
        print(f"W={window}", flush=True)
        for policy in POLICIES:
            rows = [r for r in records if r["window"] == window and r["policy"] == policy]
            times = [r["elapsed_ms"] for r in rows]
            cpu_times = [r["cpu_ms"] for r in rows]
            aggregate = {**{k: rows[0][k] for k in
                           ("window", "policy", "cp95_m", "removal_calls",
                            "max_solve_states", "final_active_states", "total_gn_iterations",
                            "max_boundary_dim", "max_local_rows")},
                         "elapsed_ms": float(np.median(times)),
                         "elapsed_q25_ms": float(np.percentile(times, 25)),
                         "elapsed_q75_ms": float(np.percentile(times, 75)),
                         "cpu_ms": float(np.median(cpu_times)),
                         "cpu_q25_ms": float(np.percentile(cpu_times, 25)),
                         "cpu_q75_ms": float(np.percentile(cpu_times, 75)),
                         "timing_suspect_runs": sum(r["timing_suspect"] for r in rows)}
            if args.profile_stages:
                for stage in ("add_state_factor_ms", "estimate_ms", "marginalize_ms", "other_ms"):
                    aggregate[stage.replace("_ms", "_mean_ms")] = float(np.mean([r[stage] for r in rows]))
                aggregate["elapsed_mean_ms"] = float(np.mean(times))
            aggregates.append(aggregate)
            print(f"  {policy:7s}: CP95={aggregate['cp95_m']:.6g} m, "
                  f"wall={aggregate['elapsed_ms']:.2f} ms, cpu={aggregate['cpu_ms']:.2f} ms, "
                  f"removals={aggregate['removal_calls']}, suspect={aggregate['timing_suspect_runs']}",
                  flush=True)
            if args.profile_stages:
                print("    stage means: " + ", ".join(
                    f"{stage}={aggregate[stage + '_mean_ms']:.2f} ms"
                    for stage in ("add_state_factor", "estimate", "marginalize", "other")), flush=True)
    if source_hashes(root) != before:
        raise RuntimeError("Experiment source changed during experiment")
    suspect = [{k: row[k] for k in ("repeat", "window", "policy", "elapsed_ms", "cpu_ms", "started_utc")}
               for row in records if row["timing_suspect"]]
    batch_control = None
    if 100 in windows:
        np.testing.assert_array_equal(trajectories["schur_w100"], trajectories["discard_w100"])
        batch_rows = {r["policy"]: r for r in aggregates if r["window"] == 100}
        batch_control = {"trajectories_identical": True,
                         "removal_calls": {p: r["removal_calls"] for p, r in batch_rows.items()},
                         "schur_vs_discard_wall_percent": 100 * (batch_rows["schur"]["elapsed_ms"] / batch_rows["discard"]["elapsed_ms"] - 1),
                         "schur_vs_discard_cpu_percent": 100 * (batch_rows["schur"]["cpu_ms"] / batch_rows["discard"]["cpu_ms"] - 1)}
    summary = {"created_utc": datetime.now(timezone.utc).isoformat(),
               "source_commit": subprocess.check_output(
                   ["git", "rev-parse", "HEAD"], cwd=root, text=True).strip(),
               "windows": sorted(args.windows), "repeats": args.repeats, "warmups": args.warmups,
               "profile_stages": args.profile_stages,
               "epochs": 100, "seed": 7, "blas_threads": 1, "python": sys.version,
               "logical_cpu": args.cpu,
               "run_order": "Windows shuffled within each repeat (seed 730); paired policies alternate first position across repeats",
               "platform": platform.platform(), "numpy": np.__version__,
               "initial_error": [100, -100, 0, 0], "max_iteration": 10, "robust_kernel": "none",
               "imitate_kfv": False,
               "data": {"anchor_radius": 105, "outlier_weight": 0., "outlier_sigma": 10., "white_sigma": .1},
               "marginalization_algorithm": "Incident factors only; pivoted local QR elimination and boundary-only square-root prior",
               "trajectory_semantics": "State at retirement; remaining window at final solve; all 100 epochs",
               "timing_scope": "Complete estimator run; excludes input generation, metrics and plotting",
               "cpu_timing_scope": "Process user plus system CPU time over the same estimator call; excludes time not executing",
               "timing_quality": {"rule": "Flag when wall minus CPU time exceeds 1000 ms AND wall time exceeds twice CPU time",
                                  "suspect_runs": suspect, "suspect_run_count": len(suspect),
                                  "measurements_discarded": 0,
                                  "batch_control": batch_control},
               "stage_timing_scope": {
                   "add_state_factor_ms": "Initial graph/prior construction, motion prediction, state/factor construction and insertion",
                   "estimate_ms": "Entire graph.estimate call including linearization, assembly, solve and state updates",
                   "marginalize_ms": "Entire marginalization/discard call including any internal assembly, prior creation and retirement; no double counting",
                   "other_ms": "Unassigned loop, output and profiling overhead; wall time minus the three stages",
                   "pie_denominator": "Sum of mean times of the three measured phases; excludes other_ms"},
               "window_behavior": "Remove oldest state only on overflow; W=100 performs zero removals",
               "retained_behavior": ["Optimization before window cleanup (may solve W+1 states)",
                   "Original motion model, noise covariance, observation schedule and final-trajectory metric"],
               "shared_fixes": ["Position prior uses positive identity Jacobian for A*delta=b, x+=delta",
                   "No forced initial-state removal; both policies release retired factors"],
               "source_hashes": before, "aggregates": aggregates}
    write_csv(output / "window_summary.csv", aggregates)
    np.savez_compressed(output / "trajectories.npz", truth=data["true_positions"], **trajectories)
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    if args.plot:
        render_plots(output, summary)
    if suspect:
        print(f"TIMING QUALITY ALERT: {len(suspect)} possible interrupted runs; all samples retained for inspection.", flush=True)
    print(f"OUTPUT: {output}", flush=True)


if __name__ == "__main__":
    main()
