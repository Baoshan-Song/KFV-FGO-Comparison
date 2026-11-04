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
def test_simulation_solvers_are_independent(mode, monkeypatch):
    cfg = load_config(resource_path("kfv_fgo_comparison.json"))
    cfg.kfv.mode = mode
    cfg = convert_kfv_to_fgo(cfg)
    data = load_data(ROOT / "data" / cfg.data_path)
    k = KfvEstimator(cfg, data).run()
    monkeypatch.setattr(KfvEstimator, "run", lambda *a: (_ for _ in ()).throw(AssertionError("delegated")))
    f = FgoEstimator(cfg, data).run()
    np.testing.assert_allclose(k["X_online"], f["X_online"], atol=1e-8, rtol=0)


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


def test_marginalization_keeps_anchor_and_affine_term():
    state = State(1, 1, np.zeros(4))
    prior = MarginFactor([state], np.eye(4), np.arange(4), np.zeros(4))
    state.value += 100
    np.testing.assert_array_equal(prior.evaluate().b, np.arange(4)-100)


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
def test_real_independent_qr_equivalence(real_data, mode, monkeypatch):
    cfg = real_config()
    cfg.kfv.mode = mode
    cfg = convert_kfv_to_fgo(cfg)
    k = KfvEstimator(cfg, real_data).run()
    monkeypatch.setattr(KfvEstimator, "run", lambda *a: (_ for _ in ()).throw(AssertionError("delegated")))
    f = FgoEstimator(cfg, real_data).run()
    assert f["solver"] == "FgoEstimator"
    np.testing.assert_allclose(k["X_online"], f["X_online"], rtol=0, atol=2e-6)
    np.testing.assert_allclose(k["covariances"], f["covariances"], rtol=1e-7, atol=2e-6)
    # EKF has no iterated-prior ambiguity and remains a useful legacy check.
    if mode == "EKF":
        np.testing.assert_allclose(k["X"], REFERENCE["real_EKF_kfv"], rtol=0, atol=1e-6)
    for result in (k,f):
        timings = result['epoch_runtime_ms']
        assert len(timings) == real_data.num_steps
        assert np.isnan(timings[0])
        assert np.all(np.isfinite(timings[1:]) & (timings[1:] > 0))
        assert np.sum(timings[1:]) <= result['runtime_seconds']*1000
        for i in (1,100,437):
            residual,_,_,ids = real_data.linearize(result["X_online"][:,i],i)
            record=result["epoch_diagnostics"][i]
            np.testing.assert_array_equal(record["postfit_ids"],ids)
            np.testing.assert_allclose(record["iteration_trace"][-1]["measurement_residual_norm_m"],np.linalg.norm(residual))


@pytest.mark.parametrize("window", [1,5])
def test_real_sliding_graph_preserves_causal_output(real_data,window,monkeypatch):
    cfg=real_config();cfg.kfv.mode="riEKF";cfg=convert_kfv_to_fgo(cfg)
    cfg.fgo.imitate_kfv=False;cfg.fgo.window_size=window
    monkeypatch.setattr(KfvEstimator,"run",lambda *a: (_ for _ in ()).throw(AssertionError("delegated")))
    result=FgoEstimator(cfg,real_data).run()
    assert result["X_online"].shape==(10,438)
    assert np.isfinite(result["X_online"]).all()
    for i,d in enumerate(result["epoch_diagnostics"][1:],1):
        assert d["active_states"]<=window+1
        np.testing.assert_array_equal(d["iteration_trace"][-1]["state"],result["X_online"][:,i])


def test_missing_files_fail_explicitly(tmp_path):
    from kfv_fgo.data.gnss_ins import GnssImuDataset

    cfg = real_config()
    cfg.data_path = str(tmp_path)
    with pytest.raises(FileNotFoundError, match="f9p_navi.obs"):
        GnssImuDataset(cfg)


def test_article_groups_isolate_loss_and_window():
    from real_experiment import article_methods
    methods = article_methods()
    assert len(methods) == 6
    for group,mode,kernel in [(methods[:3],'iEKF','none'),(methods[3:],'riEKF','huber')]:
        k,re,sw=[m['config'] for m in group]
        assert k.kfv.mode==mode
        assert re.fgo.imitate_kfv and re.fgo.window_size==1
        assert not sw.fgo.imitate_kfv and sw.fgo.window_size==20
        for cfg in (k.kfv,re.fgo,sw.fgo):
            assert cfg.max_iteration==5 and cfg.robust_kernel==kernel
            np.testing.assert_array_equal(cfg.q,k.kfv.q)
            np.testing.assert_array_equal(cfg.p0,k.kfv.p0)


def test_compact_imu_preserves_csv_values(tmp_path):
    from kfv_fgo.data.gnss_ins import ImuData
    values=np.zeros((5,32))
    values[:,0]=np.array([1,2,2,3,4])*1e9
    values[:,7]=1
    values[:,29:32]=np.arange(15).reshape(5,3)/10
    csv=tmp_path/'imu.csv';npz=tmp_path/'imu.npz'
    np.savetxt(csv,values,delimiter=',',header='columns',comments='')
    np.savez_compressed(npz,time_s=values[:,0]*1e-9,quat_xyzw=values[:,4:8],acc_mps2=values[:,29:32])
    a,b=ImuData(csv),ImuData(npz)
    for key in ('time','quat','acc'):np.testing.assert_array_equal(getattr(a,key),getattr(b,key))
    for key in ('time','quat','acc'):np.testing.assert_array_equal(a.interval(1.5,3.5)[key],b.interval(1.5,3.5)[key])


def test_real_viewer_excludes_initial_prior(real_data):
    from copy import copy
    from real_experiment import preset_pair,run_real_pair,real_payload
    short=copy(real_data);short.num_steps=4;short.timestamps=real_data.timestamps[:4]
    result=run_real_pair(short,*preset_pair('Matched EKF'))
    payload=real_payload(result)
    for item,err in zip(result['estimators'],payload['errors']):
        assert not item['errors']['valid'][0]
        assert item['accuracy']['initial_prior_excluded']
        assert item['accuracy']['evaluated_epochs']==3
        assert err['horizontal'][0] is None


def test_insufficient_satellites_has_epoch_context(real_data):
    from kfv_fgo.model.gnss_ins import GnssObservationModel

    model = GnssObservationModel(real_data.backend, {"minimum_snr": 1000})
    with pytest.raises(ValueError, match="epoch 2.*found 0"):
        model.linearize(real_data.initial_state, real_data.epochs[1])
