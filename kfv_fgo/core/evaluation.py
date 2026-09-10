"""Timestamp-aligned accuracy and algorithm-difference reports."""

import json
from pathlib import Path

import numpy as np

from ..model.gnss_ins import local_to_ecef


def statistics(error):
    error = np.asarray(error, dtype=float)
    if error.size == 0 or not np.all(np.isfinite(error)):
        raise ValueError("No finite errors to evaluate")
    return {
        "count": int(error.size),
        "mse": float(np.mean(error**2)),
        "rmse": float(np.sqrt(np.mean(error**2))),
        "mae": float(np.mean(np.abs(error))),
        "max_error": float(np.max(np.abs(error))),
        "absolute_error_95": float(np.percentile(np.abs(error), 95)),
    }


def load_ground_truth(path):
    raw = np.loadtxt(path, skiprows=2)
    if raw.ndim != 2 or raw.shape[1] < 10:
        raise ValueError("Unsupported gt.txt format")

    def dms(parts):
        sign = np.where(parts[:, 0] < 0, -1, 1)
        return sign * (np.abs(parts[:, 0]) + parts[:, 1] / 60 + parts[:, 2] / 3600)

    lat, lon = np.deg2rad(dms(raw[:, 3:6])), np.deg2rad(dms(raw[:, 6:9]))
    e2 = (1 / 298.257223563) * (2 - 1 / 298.257223563)
    n = 6378137.0 / np.sqrt(1 - e2 * np.sin(lat) ** 2)
    height = raw[:, 9]
    xyz = np.column_stack(
        (
            (n + height) * np.cos(lat) * np.cos(lon),
            (n + height) * np.cos(lat) * np.sin(lon),
            (n * (1 - e2) + height) * np.sin(lat),
        )
    )
    times = raw[:, 0]
    if (
        not np.all(np.isfinite(xyz))
        or not np.all(np.isfinite(times))
        or np.any(np.diff(times) <= 0)
    ):
        raise ValueError("Ground truth must be finite and have increasing timestamps")
    return times, xyz


def evaluate_ground_truth(result, truth):
    times, xyz = truth
    query = np.asarray(result["timestamps"])
    if result["X"].shape[1] != len(query):
        raise ValueError("State/timestamp lengths differ")
    valid = (query >= times[0]) & (query <= times[-1])
    selected = query[valid]
    if not len(selected):
        raise ValueError("No overlapping GNSS and ground-truth timestamps")
    reference = np.column_stack(
        [np.interp(selected, times, xyz[:, i]) for i in range(3)]
    )
    error_xyz = result["X"][:3, valid].T - reference
    error_enu = error_xyz @ local_to_ecef(reference[0])
    horizontal = np.linalg.norm(error_enu[:, :2], axis=1)
    spatial = np.linalg.norm(error_enu, axis=1)
    report = {
        "horizontal": statistics(horizontal),
        "three_dimensional": statistics(spatial),
        "ENU": {
            name: statistics(error_enu[:, i])
            for i, name in enumerate(("east", "north", "up"))
        },
        "evaluated_epochs": int(valid.sum()),
        "excluded_epochs": int((~valid).sum()),
        "alignment": "Linear interpolation in ECEF on overlapping UTC times; no extrapolation",
    }
    return report, {
        "timestamps": selected,
        "error_enu": error_enu,
        "horizontal": horizontal,
        "three_dimensional": spatial,
        "truth_ecef": reference,
        "valid": valid,
    }


def compare_results(left, right):
    lt, rt = np.asarray(left["timestamps"]), np.asarray(right["timestamps"])
    if left["X"].shape != right["X"].shape or not np.array_equal(lt, rt):
        raise ValueError(
            "Algorithm comparison requires identical state dimensions and timestamps"
        )
    difference = left["X"] - right["X"]
    position_dim = 3 if difference.shape[0] == 10 else 2
    return {
        "max_absolute_state_difference": float(np.max(np.abs(difference))),
        "per_state_max_absolute_difference": np.max(
            np.abs(difference), axis=1
        ).tolist(),
        "position_difference_m": statistics(
            np.linalg.norm(difference[:position_dim], axis=0)
        ),
    }


def save_results(output, results, truth=None, simulation_truth=None, plots=True):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    summary = {"algorithms": {}, "comparisons": {}}
    evaluations = {}
    for name, result in results.items():
        np.savez_compressed(
            output / f"{name}.npz", X=result["X"], timestamps=result["timestamps"]
        )
        np.savetxt(
            output / f"{name}.csv",
            np.column_stack((result["timestamps"], result["X"].T)),
            delimiter=",",
            header="utc_or_sim_seconds,"
            + ",".join(f"state_{i}" for i in range(result["X"].shape[0])),
            comments="",
        )
        report = {
            "runtime_seconds": result["runtime_seconds"],
            "solver": result.get("solver", "unknown"),
            "config": result["config"],
            "backend": result["backend"],
            "epochs": len(result["timestamps"]),
            "trajectory_semantics": result.get(
                "trajectory_semantics", "Causal filter states"
            ),
        }
        if truth is not None:
            report["accuracy"], evaluations[name] = evaluate_ground_truth(result, truth)
            e = evaluations[name]
            np.savetxt(
                output / f"{name}_errors.csv",
                np.column_stack(
                    (
                        e["timestamps"],
                        e["error_enu"],
                        e["horizontal"],
                        e["three_dimensional"],
                    )
                ),
                delimiter=",",
                header="utc_seconds,east_m,north_m,up_m,horizontal_m,three_dimensional_m",
                comments="",
            )
        if simulation_truth is not None:
            if result["X"].shape[1] != simulation_truth.shape[1]:
                raise ValueError("Simulation truth length mismatch")
            report["accuracy"] = statistics(
                np.linalg.norm(result["X"][:2] - simulation_truth[:2], axis=0)
            )
        summary["algorithms"][name] = report
    for name in results:
        if name.startswith("KFV_"):
            peer = "ReFGO_" + name[4:]
            if peer in results:
                summary["comparisons"][name + "_vs_" + peer] = compare_results(
                    results[name], results[peer]
                )
                np.savetxt(
                    output / f"{name}_vs_{peer}.csv",
                    np.column_stack(
                        (
                            results[name]["timestamps"],
                            (results[name]["X"] - results[peer]["X"]).T,
                        )
                    ),
                    delimiter=",",
                    header="timestamp,"
                    + ",".join(
                        f"delta_state_{i}" for i in range(results[name]["X"].shape[0])
                    ),
                    comments="",
                )
    if plots:
        plot_results(output, results, evaluations, simulation_truth)
    (output / "summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    return summary


def plot_results(output, results, evaluations, simulation_truth):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(9, 6))
    origin = next(iter(results.values()))["X"][:3, 0]
    is_real = next(iter(results.values()))["config"]["data"]["mode"] == "real"
    rotation = local_to_ecef(origin) if is_real else None
    if evaluations:
        e = next(iter(evaluations.values()))
        gt = (e["truth_ecef"] - origin) @ rotation
        ax.plot(gt[:, 0], gt[:, 1], "k--", label="Ground truth", linewidth=2)
    elif simulation_truth is not None:
        ax.plot(*simulation_truth[:2], "k--", label="Ground truth")
    for name, result in results.items():
        xy = (
            (result["X"][:3].T - origin) @ rotation
            if rotation is not None
            else result["X"][:2].T
        )
        ax.plot(xy[:, 0], xy[:, 1], label=name, linewidth=1)
    ax.set(
        xlabel="East (m)" if is_real else "X (m)",
        ylabel="North (m)" if is_real else "Y (m)",
        title="Trajectory comparison",
    )
    ax.axis("equal")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(output / "trajectories.png", dpi=160)
    plt.close(fig)
    if evaluations:
        fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
        t0 = min(e["timestamps"][0] for e in evaluations.values())
        for name, e in evaluations.items():
            axes[0].plot(e["timestamps"] - t0, e["horizontal"], label=name)
            axes[1].plot(e["timestamps"] - t0, e["three_dimensional"], label=name)
        axes[0].set(
            ylabel="Horizontal error (m)", title="Accuracy against ground truth"
        )
        axes[1].set(ylabel="3D error (m)", xlabel="Time from evaluation start (s)")
        for ax in axes:
            ax.grid(alpha=0.3)
        axes[0].legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(output / "position_errors.png", dpi=160)
        plt.close(fig)
    pairs = [
        (name, "ReFGO_" + name[4:])
        for name in results
        if name.startswith("KFV_") and "ReFGO_" + name[4:] in results
    ]
    if pairs:
        fig, ax = plt.subplots(figsize=(10, 4))
        for name, peer in pairs:
            n = 3 if is_real else 2
            difference = np.linalg.norm(
                results[name]["X"][:n] - results[peer]["X"][:n], axis=0
            )
            t = results[name]["timestamps"]
            ax.plot(t - t[0], difference, label=name[4:])
        ax.set(
            xlabel="Time (s)",
            ylabel="Position difference (m)",
            title="KFV versus FGO template (original KFV shortcut)",
        )
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        ax.grid(alpha=0.3)
        ax.legend()
        fig.tight_layout()
        fig.savefig(output / "algorithm_differences.png", dpi=160)
        plt.close(fig)
