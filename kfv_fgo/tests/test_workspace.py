import numpy as np
from kfv_fgo.config.settings import (
    load_config,
    resolve_data_path,
    resource_path,
    workspace_root,
)
from kfv_fgo.core.estimator import FgoEstimator, KfvEstimator
from kfv_fgo.data.simulation import load_data


def test_assets_resolve_inside_workspace():
    for filename in ("swfgo.json", "gnss_ins.json"):
        config_path = resource_path(filename)
        cfg = load_config(config_path)
        assert config_path.is_relative_to(workspace_root())
        assert resolve_data_path(cfg).is_relative_to(workspace_root())
        assert resolve_data_path(cfg).exists()


def test_sliding_graph_does_not_delegate_to_kfv(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("SWFGO must solve the factor graph")

    monkeypatch.setattr(KfvEstimator, "run", forbidden)
    cfg = load_config(resource_path("swfgo.json"))
    cfg.fgo.window_size = 5
    data = load_data(resolve_data_path(cfg))
    data["num_steps"] = 12
    for key in ("true_positions", "true_velocities", "toa_measurements"):
        data[key] = data[key][:, :12]
    result = FgoEstimator(cfg, data).run()
    assert result["solver"] == "FgoEstimator"
    assert result["X"].shape == (4, 12)
    assert np.isfinite(result["X"]).all()
