import hashlib
import json
from pathlib import Path

import numpy as np
import pytest
from kfv_fgo.examples.compare_matlab_results import compare_results, load_result, run
from scipy.io import savemat


def saved_pair(tmp_path, n=4):
    x = np.arange(n * 3, dtype=float).reshape(n, 3)
    times = np.array([0.0, 1.0, 2.0])
    matlab = tmp_path / "matlab.mat"
    python = tmp_path / "python.npz"
    savemat(matlab, {"reference": {"runs": {"KFV_EKF": {"X": x, "timestamps": times}}}})
    np.savez(python, X=x, timestamps=times)
    return matlab, python


@pytest.mark.parametrize("n", [4, 10])
def test_saved_results_compare_without_launching_any_process(tmp_path, monkeypatch, n):
    import subprocess

    def forbidden(*args, **kwargs):
        raise AssertionError(
            "A saved-result comparison must not launch MATLAB or Python estimators"
        )

    monkeypatch.setattr(subprocess, "Popen", forbidden)
    matlab, python = saved_pair(tmp_path, n)
    output = tmp_path / "comparison"
    report = run(
        matlab, python, output, matlab_key="reference.runs.KFV_EKF", plots=False
    )
    assert report["passed"]
    assert report["epochs"] == 3
    assert report["max_absolute_state_difference"] == 0
    assert json.loads((output / "metrics.json").read_text())["state_dimension"] == n
    assert np.load(output / "comparison.npz")["difference"].shape == (n, 3)


def test_original_matlab_field_and_simulation_dt(tmp_path):
    path = tmp_path / "original.mat"
    savemat(
        path,
        {
            "result_ekf": {"X": np.ones((4, 3))},
            "result_fg_ekf": {"X": np.zeros((4, 3))},
            "data": {"dt": 0.5},
        },
    )
    with pytest.raises(ValueError, match="Choose a result field.*result_ekf"):
        load_result(path)
    result = load_result(path, "result_ekf")
    np.testing.assert_array_equal(result["timestamps"], [0, 0.5, 1])
    assert result["source"]["time_source"] == "data.dt"


def test_real_result_requires_actual_timestamps(tmp_path):
    path = tmp_path / "real.mat"
    savemat(path, {"result": {"X": np.ones((10, 3))}, "gps": {"utc": [10, 11, 12]}})
    with pytest.raises(ValueError, match="require timestamps"):
        load_result(path, "result", dt=1)
    np.testing.assert_array_equal(
        load_result(path, "result", time_key="gps.utc")["timestamps"], [10, 11, 12]
    )


def test_missing_simulation_times_require_explicit_dt(tmp_path):
    path = tmp_path / "result.npz"
    np.savez(path, X=np.ones((4, 3)))
    with pytest.raises(ValueError, match="specify --dt"):
        load_result(path)
    np.testing.assert_array_equal(
        load_result(path, dt=2, start_time=5)["timestamps"], [5, 7, 9]
    )


def test_equal_positions_do_not_hide_velocity_difference(tmp_path):
    matlab, python = saved_pair(tmp_path)
    left, right = load_result(matlab), load_result(python)
    right["X"][3, 1] += 0.1
    report, _ = compare_results(left, right)
    assert not report["passed"]
    assert report["position_difference_m"]["max_absolute"] == 0
    assert report["velocity_difference_m_per_s"]["max_absolute"] > 0


def test_does_not_truncate_history_or_fit_timestamps(tmp_path):
    matlab, python = saved_pair(tmp_path)
    left, right = load_result(matlab), load_result(python)
    right["timestamps"] += 0.01
    with pytest.raises(ValueError, match="timestamps differ"):
        compare_results(left, right)
    right["X"] = right["X"][:, 1:]
    with pytest.raises(ValueError, match="histories will not be truncated"):
        compare_results(left, right)


@pytest.mark.parametrize("bad", [np.nan, np.inf])
def test_nonfinite_state_is_not_a_valid_comparison(tmp_path, bad):
    path = tmp_path / "bad.npz"
    x = np.ones((4, 3))
    x[0, 1] = bad
    np.savez(path, X=x, timestamps=[0, 1, 2])
    with pytest.raises(ValueError, match="finite"):
        load_result(path)


@pytest.mark.parametrize("mode,dimension", [("simulation", 4), ("real", 10)])
def test_published_matlab_fixtures_match_provenance_and_default_data(mode, dimension):
    root = Path(__file__).resolve().parents[1]
    fixtures = root / "tests/fixtures"
    metadata = json.loads((fixtures / "matlab_provenance.json").read_text())[mode]
    path = fixtures / metadata["file"]
    assert hashlib.sha256(path.read_bytes()).hexdigest() == metadata["sha256"]
    for name, digest in metadata["data_sha256"].items():
        assert (
            hashlib.sha256(
                (root / metadata["data_directory"] / name).read_bytes()
            ).hexdigest()
            == digest
        )
    for name, settings in metadata["algorithms"].items():
        result = load_result(path, f"reference.runs.{name}")
        assert result["X"].shape == (dimension, settings["epochs"])
