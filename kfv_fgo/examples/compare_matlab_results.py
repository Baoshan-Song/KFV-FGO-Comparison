"""Compare saved MATLAB .mat and Python .npz results; run neither estimator.

Example:
  python examples/compare_matlab_results.py --matlab matlab.mat \
      --matlab-key result_ekf --python results/KFV_EKF.npz

Both four-state simulation and ten-state GNSS/INS are supported. For simulation,
missing timestamps can be derived from data.dt/settings.dt or explicit --dt.
Real results must contain numeric timestamps in the same time scale.
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.io import loadmat


def field(value, key):
    for part in key.split("."):
        value = value[int(part)] if isinstance(value, (list, tuple)) else value[part]
    return value


def result_fields(value, prefix=""):
    found = {}
    if isinstance(value, dict):
        if "X" in value and np.size(value["X"]):
            return {prefix: value}
        for key, child in value.items():
            if not key.startswith("__"):
                found.update(result_fields(child, f"{prefix}.{key}".strip(".")))
    elif isinstance(value, (list, tuple)):
        for index, child in enumerate(value):
            found.update(result_fields(child, f"{prefix}.{index}".strip(".")))
    return found


def load_result(path, key=None, *, dt=None, start_time=0.0, time_key=None):
    """Read numeric result structs; ignore saved MATLAB estimator objects."""
    path = Path(path).expanduser().resolve()
    if path.suffix.lower() == ".npz":
        with np.load(path, allow_pickle=False) as archive:
            values = {name: archive[name] for name in archive.files}
    elif path.suffix.lower() == ".mat":
        try:
            values = loadmat(path, simplify_cells=True)
        except NotImplementedError as exc:
            raise ValueError(
                "Use a standard MAT (-v7) file; MAT v7.3 is not supported"
            ) from exc
    else:
        raise ValueError("Result files must be .mat or .npz")

    if key:
        try:
            result = field(values, key)
        except (KeyError, IndexError, TypeError, ValueError) as exc:
            raise ValueError(f"Result field {key!r} not found in {path.name}") from exc
    else:
        candidates = result_fields(values)
        if len(candidates) != 1:
            raise ValueError(
                f"Choose a result field in {path.name}; available: "
                + ", ".join(candidates)
            )
        key, result = next(iter(candidates.items()))
    if not isinstance(result, dict) or "X" not in result:
        raise ValueError("Select a result struct containing X, not the X field itself")
    x = np.asarray(result["X"], dtype=float)
    if x.ndim == 1 and x.size in (4, 10):
        x = x[:, None]
    if x.ndim != 2 or x.shape[0] not in (4, 10) or x.shape[1] == 0:
        raise ValueError(
            "X must have shape 4 x epochs or 10 x epochs; no automatic transpose"
        )

    times = result.get("timestamps")
    time_source = "timestamps"
    if time_key:
        try:
            times = field(values, time_key)
        except (KeyError, IndexError, TypeError, ValueError) as exc:
            raise ValueError(f"Time field {time_key!r} not found") from exc
        time_source = time_key
    if times is None:
        if x.shape[0] == 10:
            raise ValueError(
                "Real GNSS/INS results require timestamps; use --matlab-time-key if needed"
            )
        time_source = "explicit dt"
        if dt is None:
            for source, name in (
                (result, "settings.dt"),
                (values, "data.dt"),
                (values, "config.KFV.dt"),
            ):
                try:
                    dt = float(field(source, name))
                    time_source = name
                    break
                except (KeyError, TypeError, ValueError):
                    continue
        if dt is None or not np.isfinite(dt) or dt <= 0:
            raise ValueError(
                "Simulation result lacks timestamps and a positive dt; specify --dt"
            )
        times = start_time + np.arange(x.shape[1]) * dt
    times = np.asarray(times, dtype=float).reshape(-1)
    if times.size != x.shape[1] or not np.all(np.diff(times) > 0):
        raise ValueError(
            "Timestamps must match the number of epochs and be strictly increasing"
        )
    if not np.isfinite(x).all() or not np.isfinite(times).all():
        raise ValueError("States and timestamps must be finite")
    return {
        "X": x,
        "timestamps": times,
        "source": {
            "file": str(path),
            "field": key or "X",
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "time_source": time_source,
        },
    }


def statistics(values):
    values = np.asarray(values)
    return {
        "rmse": float(np.sqrt(np.mean(values**2))),
        "mae": float(np.mean(np.abs(values))),
        "p95": float(np.percentile(np.abs(values), 95)),
        "max_absolute": float(np.max(np.abs(values))),
    }


def compare_results(matlab, python, atol=None):
    mx, px = matlab["X"], python["X"]
    if mx.shape != px.shape:
        raise ValueError(
            "MATLAB/Python state shapes differ; histories will not be truncated"
        )
    if not np.array_equal(matlab["timestamps"], python["timestamps"]):
        raise ValueError("MATLAB/Python timestamps differ; no automatic time fitting")
    real = mx.shape[0] == 10
    if atol is None:
        atol = 1e-6 if real else 1e-8
    if not np.isfinite(atol) or atol < 0:
        raise ValueError("Tolerance must be finite and nonnegative")
    names = (
        (
            "ecef_x",
            "ecef_y",
            "ecef_z",
            "vx",
            "vy",
            "vz",
            "bias_x",
            "bias_y",
            "bias_z",
            "clock",
        )
        if real
        else ("x", "y", "vx", "vy")
    )
    difference = px - mx
    position_dim = 3 if real else 2
    report = {
        "passed": bool(np.all(np.abs(difference) <= atol)),
        "epochs": mx.shape[1],
        "state_dimension": mx.shape[0],
        "absolute_tolerance": atol,
        "relative_tolerance": 0,
        "difference_convention": "Python minus MATLAB",
        "scope": "Saved state histories only; matching model/configuration is not inferred",
        "state_units": ["m"] * position_dim
        + ["m/s"] * position_dim
        + (["m/s^2"] * 3 + ["m"] if real else []),
        "matlab": matlab["source"],
        "python": python["source"],
        "max_absolute_state_difference": float(np.max(np.abs(difference))),
        "per_state": {
            name: statistics(delta) for name, delta in zip(names, difference)
        },
        "position_difference_m": statistics(
            np.linalg.norm(difference[:position_dim], axis=0)
        ),
        "velocity_difference_m_per_s": statistics(
            np.linalg.norm(difference[position_dim : 2 * position_dim], axis=0)
        ),
    }
    return report, difference


def plot_results(output, matlab, python, report, difference):
    import matplotlib.pyplot as plt

    names = list(report["per_state"])
    elapsed = matlab["timestamps"] - matlab["timestamps"][0]
    fig, axes = plt.subplots(len(names), 1, sharex=True, figsize=(10, 2 * len(names)))
    for axis, name, delta, unit in zip(axes, names, difference, report["state_units"]):
        axis.plot(elapsed, delta)
        axis.set_ylabel(f"{name} ({unit})")
        axis.grid(True)
    axes[0].set_title("Python minus MATLAB: complete state histories")
    axes[-1].set_xlabel("Elapsed time (s)")
    fig.tight_layout()
    fig.savefig(output / "state_differences.png", dpi=150)
    plt.close(fig)
    fig, axis = plt.subplots(figsize=(7, 6))
    real = len(names) == 10
    origin = matlab["X"][:2, :1] if real else np.zeros((2, 1))
    for label, result, style in (("MATLAB", matlab, "-"), ("Python", python, "--")):
        position = result["X"][:2] - origin
        axis.plot(*position, style, label=label)
    axis.set_xlabel("ECEF X offset (m)" if real else "X (m)")
    axis.set_ylabel("ECEF Y offset (m)" if real else "Y (m)")
    axis.axis("equal")
    axis.grid(True)
    axis.legend()
    fig.tight_layout()
    fig.savefig(output / "trajectories.png", dpi=150)
    plt.close(fig)


def run(
    matlab_path,
    python_path,
    output=None,
    *,
    matlab_key=None,
    python_key=None,
    matlab_time_key=None,
    python_time_key=None,
    dt=None,
    start_time=0.0,
    atol=None,
    plots=True,
):
    matlab = load_result(
        matlab_path, matlab_key, dt=dt, start_time=start_time, time_key=matlab_time_key
    )
    python = load_result(
        python_path, python_key, dt=dt, start_time=start_time, time_key=python_time_key
    )
    report, difference = compare_results(matlab, python, atol)
    output = (
        Path(output)
        if output is not None
        else Path.cwd() / "results" / "matlab_comparison"
    )
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(
            f"Output directory is not empty: {output}; choose another --output"
        )
    output.mkdir(parents=True, exist_ok=True)
    names = list(report["per_state"])
    np.savez_compressed(
        output / "comparison.npz",
        timestamps=matlab["timestamps"],
        matlab_X=matlab["X"],
        python_X=python["X"],
        difference=difference,
    )
    header = ["timestamp"] + [
        f"{prefix}_{name}"
        for prefix in ("matlab", "python", "difference")
        for name in names
    ]
    np.savetxt(
        output / "comparison.csv",
        np.column_stack(
            (matlab["timestamps"], matlab["X"].T, python["X"].T, difference.T)
        ),
        delimiter=",",
        header=",".join(header),
        comments="",
    )
    (output / "metrics.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    if plots:
        plot_results(output, matlab, python, report, difference)
    print(
        f"{'PASS' if report['passed'] else 'DIFFERENT'}: {report['epochs']} epochs; "
        f"maximum state difference {report['max_absolute_state_difference']:.9g}"
    )
    print(f"Saved comparison to {output.resolve()}")
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--matlab",
        required=True,
        type=Path,
        help="Saved MAT file, not a MATLAB executable",
    )
    parser.add_argument(
        "--python", required=True, type=Path, help="Saved Python NPZ or MAT result"
    )
    parser.add_argument(
        "--matlab-key", help="For example result_ekf or reference.runs.KFV_EKF"
    )
    parser.add_argument("--python-key")
    parser.add_argument(
        "--matlab-time-key", help="Numeric timestamp field, including dotted paths"
    )
    parser.add_argument("--python-time-key")
    parser.add_argument(
        "--dt", type=float, help="Used only for simulation results without timestamps"
    )
    parser.add_argument("--start-time", type=float, default=0.0)
    parser.add_argument("--atol", type=float)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--no-plots", action="store_true")
    args = parser.parse_args(argv)
    try:
        report = run(
            args.matlab,
            args.python,
            args.output,
            matlab_key=args.matlab_key,
            python_key=args.python_key,
            matlab_time_key=args.matlab_time_key,
            python_time_key=args.python_time_key,
            dt=args.dt,
            start_time=args.start_time,
            atol=args.atol,
            plots=not args.no_plots,
        )
    except (ValueError, OSError) as exc:
        parser.exit(1, f"Comparison failed: {exc}\n")
    return 0 if report["passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
