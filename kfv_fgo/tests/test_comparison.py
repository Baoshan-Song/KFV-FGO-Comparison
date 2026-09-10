from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest
from kfv_fgo.config.settings import (
    convert_kfv_to_fgo,
    load_config,
    real_config,
    resource_path,
    save_config,
)
from kfv_fgo.core.estimator import FgoEstimator, KfvEstimator
from kfv_fgo.core.evaluation import compare_results, evaluate_ground_truth
from kfv_fgo.core.fgo import Factor, FactorGraph, MarginFactor, State
from kfv_fgo.data.simulation import load_data
from kfv_fgo.model.gnss_ins import ImuPropagationModel

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = np.load(Path(__file__).parent / "fixtures" / "matlab_reference.npz")
MODES = ("EKF", "iEKF", "rEKF", "riEKF")


@pytest.mark.parametrize("mode", MODES)
def test_simulation_matches_independent_matlab_reference(mode, monkeypatch):
    cfg = load_config(resource_path("kfv_fgo_comparison.json"))
    cfg.kfv.mode = mode
    data = load_data(ROOT / "data" / cfg.data_path)
    k = KfvEstimator(cfg, data).run()
    np.testing.assert_allclose(k["X"], REFERENCE[f"sim_{mode}_kfv"], rtol=0, atol=1e-8)

    original_run = KfvEstimator.run
    calls = []

    def delegated(estimator):
        calls.append(estimator)
        return original_run(estimator)

    monkeypatch.setattr(KfvEstimator, "run", delegated)
    f = FgoEstimator(convert_kfv_to_fgo(cfg), data).run()
    assert len(calls) == 1
    np.testing.assert_array_equal(f["X"], k["X"])
    assert f["solver"] == "KfvEstimator"
    np.testing.assert_allclose(f["X"], REFERENCE[f"sim_{mode}_fgo"], rtol=0, atol=1e-8)
    for i in range(1, 4):
        for key, original in [
            ("residual", "residual"),
            ("kalman_gain", "Kalman_gain"),
            ("innovation_covariance", "innovation_covariance"),
        ]:
            np.testing.assert_allclose(
                k["debug_info"][i][key],
                REFERENCE[f"sim_{mode}_debug{i}_{original}"],
                rtol=1e-9,
                atol=1e-8,
            )


def test_sliding_returns_every_epoch_and_honors_window():
    cfg = load_config(resource_path("swfgo.json"))
    cfg.fgo.window_size = 3
    data = load_data(ROOT / "data" / cfg.data_path)
    data["num_steps"] = 12
    data["true_positions"] = data["true_positions"][:, :12]
    data["true_velocities"] = data["true_velocities"][:, :12]
    data["toa_measurements"] = data["toa_measurements"][:, :12]
    result = FgoEstimator(cfg, data).run()
    assert result["X"].shape == (4, 12)
    assert all(d["active_states"] <= 3 for d in result["debug_info"][1:])
    assert np.array_equal(result["timestamps"], np.arange(12))


def test_normal_equation_uses_actual_residual_sizes():
    class Fixed(Factor):
        def evaluate(self):
            return self

    graph = FactorGraph(real_config().fgo)
    state = State(1, 1, np.zeros(10))
    graph.add_state(state)
    for rows in (1, 4, 7):
        f = Fixed([state])
        f.A = np.ones((rows, 10))
        f.b = np.arange(rows)
        graph.add_factor(f)
    graph.normal_equation()
    assert graph.J.shape == (12, 10)
    assert graph.r.shape == (12,)


def test_marginalization_preserves_matlab_frozen_prior():
    state = State(1, 1, np.zeros(4))
    prior = MarginFactor([state], np.eye(4), np.arange(4))
    state.value += 100
    np.testing.assert_array_equal(prior.evaluate().b, np.arange(4))


def test_imu_keeps_first_sample_euler_convention():
    x = np.array([6378137.0, 0, 0, 1, 2, 3, 0, 0, 0, 7])
    batch = {
        "time": np.array([100, 102]),
        "acc": np.array([[0, 0, 9.81], [5, 5, 5]]),
        "quat": np.array([[0, 0, 0, 1], [0, 0, 0, 1]]),
    }
    prediction, F, Q = ImuPropagationModel(np.eye(10)).propagate(x, batch)
    np.testing.assert_allclose(prediction[:3], x[:3] + 2 * x[3:6], rtol=0, atol=1e-10)
    np.testing.assert_allclose(prediction[3:], x[3:], rtol=0, atol=1e-10)
    np.testing.assert_array_equal(F[:3, 3:6], 2 * np.eye(3))
    batch["acc"][:] = np.nan
    with pytest.raises(ValueError, match="no valid"):
        ImuPropagationModel(Q).propagate(x, batch)


def test_accuracy_interpolates_time_without_extrapolation():
    origin = np.array([6378137.0, 0, 0])
    truth = (
        np.array([10.0, 11.0, 12.0]),
        origin + np.array([[0, 0, 0], [0, 2, 0], [0, 4, 0]]),
    )
    times = np.array([9.0, 10.5, 11.5, 13.0])
    result = {
        "timestamps": times,
        "X": (origin + np.array([[0, -2, 0], [0, 1, 0], [0, 3, 0], [0, 6, 0]])).T,
    }
    report, errors = evaluate_ground_truth(result, truth)
    assert report["excluded_epochs"] == 2
    assert report["horizontal"]["rmse"] == 0
    np.testing.assert_array_equal(errors["timestamps"], [10.5, 11.5])
    other = deepcopy(result)
    other["timestamps"] = times + 0.001
    with pytest.raises(ValueError, match="timestamps"):
        compare_results(result, other)


def test_real_config_roundtrip(tmp_path):
    config = real_config()
    path = tmp_path / "real.json"
    save_config(config, path)
    loaded = load_config(path)
    assert loaded.state_dim == 10
    assert loaded.gnss == config.gnss
    assert loaded.kfv.r is None
    np.testing.assert_array_equal(loaded.kfv.q, config.kfv.q)


def test_utf8_json_path(tmp_path):
    import json

    path = tmp_path / "settings.json"
    data_path = "\u6570\u636e/\u57ce\u5e02"
    path.write_text(
        json.dumps({"data": {"path": data_path}}, ensure_ascii=False), encoding="utf-8"
    )
    assert load_config(path).data_path == data_path


@pytest.fixture(scope="module")
def real_data():
    pytest.importorskip("pyrtklib")
    from kfv_fgo.data.gnss_ins import GnssImuDataset

    return GnssImuDataset(real_config())


def test_real_model_matches_matlab_probes_and_is_state_independent(real_data):
    data = real_data
    assert data.num_steps == 438
    np.testing.assert_allclose(
        data.initial_state, REFERENCE["initial_state"], rtol=0, atol=1e-6
    )
    for i in range(1, 4):
        xp, F, Q = data.propagate(REFERENCE[f"probe{i}_x"], i, real_config().kfv)
        np.testing.assert_allclose(xp, REFERENCE[f"probe{i}_xp"], rtol=0, atol=1e-7)
        np.testing.assert_allclose(F, REFERENCE[f"probe{i}_F"], rtol=0, atol=1e-12)
        packet = data.get_measurement(i)
        snapshot = packet.positions.copy()
        r, H, R, ids = data.linearize(REFERENCE[f"probe{i}_xp"], i)
        np.testing.assert_array_equal(ids, REFERENCE[f"probe{i}_sat_ids"])
        np.testing.assert_allclose(
            r, REFERENCE[f"probe{i}_residual"], rtol=0, atol=1e-6
        )
        np.testing.assert_allclose(H, REFERENCE[f"probe{i}_H"], rtol=0, atol=1e-12)
        np.testing.assert_allclose(R, REFERENCE[f"probe{i}_R"], rtol=1e-10, atol=1e-10)
        shifted = xp.copy()
        shifted[:3] += 1000
        changed = data.linearize(shifted, i)
        assert not np.allclose(changed[0], r)
        np.testing.assert_array_equal(packet.positions, snapshot)
        repeated = data.linearize(REFERENCE[f"probe{i}_xp"], i)
        np.testing.assert_array_equal(repeated[0], r)


@pytest.mark.parametrize("mode", MODES)
def test_all_real_kfv_and_recursive_fgo_match_matlab(real_data, mode, monkeypatch):
    cfg = real_config()
    cfg.kfv.mode = mode
    k = KfvEstimator(cfg, real_data).run()

    original_run = KfvEstimator.run
    calls = []

    def delegated(estimator):
        calls.append(estimator)
        return original_run(estimator)

    monkeypatch.setattr(KfvEstimator, "run", delegated)
    f = FgoEstimator(convert_kfv_to_fgo(cfg), real_data).run()
    assert len(calls) == 1
    np.testing.assert_array_equal(f["X"], k["X"])
    for name, result in [("kfv", k), ("fgo", f)]:
        assert result["X"].shape == (10, 438)
        np.testing.assert_allclose(
            result["X"], REFERENCE[f"real_{mode}_{name}"], rtol=0, atol=1e-6
        )
    for i in range(1, 4):
        for key, original in [
            ("jacobian_all", "jacobian_all"),
            ("residual", "residual"),
            ("kalman_gain", "Kalman_gain"),
            ("innovation_covariance", "innovation_covariance"),
        ]:
            np.testing.assert_allclose(
                k["debug_info"][i][key],
                REFERENCE[f"real_{mode}_debug{i}_{original}"],
                rtol=1e-7,
                atol=1e-6,
            )


@pytest.mark.parametrize("window", [1, 5, 10])
def test_full_real_sliding_matlab_parity(real_data, window):
    cfg = real_config()
    cfg.kfv.mode = "riEKF"
    cfg = convert_kfv_to_fgo(cfg)
    cfg.fgo.imitate_kfv = False
    cfg.fgo.window_size = window
    result = FgoEstimator(cfg, real_data).run()
    assert result["X"].shape == (10, 438)
    np.testing.assert_allclose(
        result["X"], REFERENCE[f"sliding_w{window}"], rtol=0, atol=1e-6
    )


def test_missing_files_fail_explicitly(tmp_path):
    from kfv_fgo.data.gnss_ins import GnssImuDataset

    cfg = real_config()
    cfg.data_path = str(tmp_path)
    with pytest.raises(FileNotFoundError, match="f9p_navi.obs"):
        GnssImuDataset(cfg)


def test_insufficient_satellites_has_epoch_context(real_data):
    from kfv_fgo.model.gnss_ins import GnssObservationModel

    model = GnssObservationModel(real_data.backend, {"minimum_snr": 1000})
    with pytest.raises(ValueError, match="epoch 2.*found 0"):
        model.linearize(real_data.initial_state, real_data.epochs[1])


def test_legacy_numerical_robust_factor_supports_vector_residuals():
    from kfv_fgo.core.fgo import AutoDiffRobustFactor

    state = State(1, 1, np.array([1.0, 2.0]))
    factor = AutoDiffRobustFactor([state], np.zeros(2), lambda x, z: x - z, lambda r: r)
    factor.evaluate()
    np.testing.assert_allclose(factor.A, -np.eye(2), rtol=0, atol=1e-10)
    np.testing.assert_array_equal(factor.b, [1.0, 2.0])
