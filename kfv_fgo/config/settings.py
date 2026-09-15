from __future__ import annotations

import json
from copy import deepcopy
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np


@dataclass
class FilterConfig:
    mode: str = "EKF"
    dt: float = 1.0
    err_x0: np.ndarray = field(
        default_factory=lambda: np.array([100.0, -100.0, 0.0, 0.0])
    )
    p0: np.ndarray = field(default_factory=lambda: np.diag([50.0, 50.0, 1.0, 1.0]))
    q: np.ndarray = field(default_factory=lambda: np.eye(4) * 1e-4)
    omega: float = 2 * np.pi / 100
    r: float = 0.1**2
    max_iteration: int = 2
    threshold_iteration: float = 1e-6
    robust_kernel: str = "huber"
    robust_delta: float = 2.0
    window_size: int = 1


@dataclass
class FgoConfig:
    dt: float = 1.0
    err_x0: np.ndarray = field(
        default_factory=lambda: np.array([100.0, -100.0, 0.0, 0.0])
    )
    p0: np.ndarray = field(default_factory=lambda: np.diag([50.0, 50.0, 1.0, 1.0]))
    q: np.ndarray = field(default_factory=lambda: np.eye(4) * 1e-4)
    omega: float = 2 * np.pi / 100
    r: float = 0.1**2
    max_iteration: int = 1
    threshold_iteration: float = 1e-6
    robust_kernel: str = "huber"
    robust_delta: float = 2.0
    window_size: int = 1
    imitate_kfv: bool = True
    autodiff: bool = False


@dataclass
class Config:
    data_mode: str = "sim"
    data_path: str = "circle_cv_gmm_L4.mat"
    state_dim: int = 4
    gnss: dict = field(default_factory=dict)
    imu: dict = field(default_factory=dict)
    leap_seconds: int = 18
    kfv: FilterConfig = field(default_factory=FilterConfig)
    fgo: FgoConfig = field(default_factory=FgoConfig)


def _config_to_dict(config: Config) -> dict:
    return {
        "state_dim": config.state_dim,
        "GNSS": config.gnss,
        "IMU": config.imu,
        "data": {
            "mode": config.data_mode,
            "path": config.data_path,
            "leap_seconds": config.leap_seconds,
        },
        "KFV": {
            "mode": config.kfv.mode,
            "dt": config.kfv.dt,
            "errX0": config.kfv.err_x0.tolist(),
            "P0": config.kfv.p0.tolist(),
            "Q": config.kfv.q.tolist(),
            "omega": config.kfv.omega,
            "R": config.kfv.r,
            "max_iteration": config.kfv.max_iteration,
            "thres_iteration": config.kfv.threshold_iteration,
            "robust_kernel": config.kfv.robust_kernel,
            "robust_delta": config.kfv.robust_delta,
            "window_size": config.kfv.window_size,
        },
        "FGO": {
            "dt": config.fgo.dt,
            "errX0": config.fgo.err_x0.tolist(),
            "P0": config.fgo.p0.tolist(),
            "Q": config.fgo.q.tolist(),
            "omega": config.fgo.omega,
            "R": config.fgo.r,
            "max_iteration": config.fgo.max_iteration,
            "thres_iteration": config.fgo.threshold_iteration,
            "robust_kernel": config.fgo.robust_kernel,
            "robust_delta": config.fgo.robust_delta,
            "window_size": config.fgo.window_size,
            "imitate_KFV": config.fgo.imitate_kfv,
            "autoDiff": config.fgo.autodiff,
        },
    }


def save_config(config: Config, path: str | Path) -> None:
    Path(path).write_text(
        json.dumps(_config_to_dict(config), indent=2) + "\n", encoding="utf-8"
    )


def load_config(path: str | Path) -> Config:
    raw = json.loads(Path(path).read_text(encoding="utf-8"))
    data = raw.get("data", {})
    k = raw.get("KFV", {})
    f = raw.get("FGO", {})
    kfv = FilterConfig(
        mode=k.get("mode", "EKF"),
        dt=k.get("dt", 1.0),
        err_x0=np.asarray(k.get("errX0", [100, -100, 0, 0]), dtype=float),
        p0=np.asarray(k.get("P0", np.diag([50, 50, 1, 1])), dtype=float),
        q=np.asarray(k.get("Q", np.eye(4) * 1e-4), dtype=float),
        omega=k.get("omega", 2 * np.pi / 100),
        r=k.get("R", 0.01),
        max_iteration=k.get("max_iteration", 2),
        threshold_iteration=k.get("thres_iteration", 1e-6),
        robust_kernel=k.get("robust_kernel", "huber"),
        robust_delta=k.get("robust_delta", 2.0),
        window_size=k.get("window_size", 1),
    )
    fgo = FgoConfig(
        dt=f.get("dt", kfv.dt),
        err_x0=np.asarray(f.get("errX0", kfv.err_x0), dtype=float),
        p0=np.asarray(f.get("P0", kfv.p0), dtype=float),
        q=np.asarray(f.get("Q", kfv.q), dtype=float),
        omega=f.get("omega", kfv.omega),
        r=f.get("R", kfv.r),
        max_iteration=f.get("max_iteration", 1),
        threshold_iteration=f.get("thres_iteration", 1e-6),
        robust_kernel=f.get("robust_kernel", "none"),
        robust_delta=f.get("robust_delta", 2.0),
        window_size=f.get("window_size", 1),
        imitate_kfv=f.get("imitate_KFV", False),
        autodiff=f.get("autoDiff", False),
    )
    return Config(
        data_mode=data.get("mode", "sim"),
        data_path=data.get("path", "circle_cv_gmm_L4.mat"),
        kfv=kfv,
        fgo=fgo,
        state_dim=raw.get("state_dim", 4),
        gnss=raw.get("GNSS", {}),
        imu=raw.get("IMU", {}),
        leap_seconds=data.get("leap_seconds", 18),
    )


def comparison_config() -> Config:
    return Config(kfv=FilterConfig(), fgo=FgoConfig(imitate_kfv=True, autodiff=False))


def sw_fgo_config() -> Config:
    return Config(
        fgo=FgoConfig(
            imitate_kfv=False,
            autodiff=False,
            max_iteration=1,
            robust_kernel="none",
            window_size=1,
        )
    )


def convert_kfv_to_fgo(config: Config) -> Config:
    k = config.kfv
    f = FgoConfig(
        dt=k.dt,
        err_x0=k.err_x0.copy(),
        p0=k.p0.copy(),
        q=k.q.copy(),
        omega=k.omega,
        r=k.r,
        max_iteration=k.max_iteration,
        threshold_iteration=k.threshold_iteration,
        robust_kernel=k.robust_kernel,
        robust_delta=k.robust_delta,
        window_size=k.window_size,
        imitate_kfv=True,
        autodiff=False,
    )
    if k.mode == "EKF":
        f.max_iteration, f.robust_kernel, f.window_size = 1, "none", 1
    elif k.mode == "iEKF":
        f.robust_kernel, f.window_size = "none", 1
    elif k.mode == "rEKF":
        f.max_iteration, f.window_size = 1, 1
    elif k.mode == "riEKF":
        f.window_size = 1
    else:
        raise ValueError(f"Unsupported KFV mode: {k.mode}")
    result = deepcopy(config)
    result.fgo = f
    return result


def real_config() -> Config:
    """Defaults from MATLAB init_settings_gnss_ins.m; no implicit gating."""
    k = FilterConfig(
        err_x0=np.zeros(10),
        p0=np.diag(np.array([5] * 6 + [500] * 3 + [100], dtype=float) ** 2),
        q=np.diag(np.array([0.3] * 3 + [0.15] * 3 + [0.01] * 3 + [100]) ** 2),
        r=None,
        max_iteration=5,
        robust_delta=1,
    )
    cfg = Config(
        data_mode="real",
        data_path="urban_nav_deep",
        state_dim=10,
        kfv=k,
        gnss={
            "minimum_snr": None,
            "minimum_elevation_deg": None,
            "Tref": 45,
            "a": 30,
            "A": 30,
            "Fref": 10,
        },
        imu={"gravity_norm": 9.81},
    )
    return convert_kfv_to_fgo(cfg)


def validate_config(config, initial, cfg):
    n = len(initial)
    if n != config.state_dim:
        raise ValueError(
            f"Initial state has {n} elements; state_dim={config.state_dim}"
        )
    for name in ("p0", "q"):
        matrix = np.asarray(getattr(cfg, name))
        if matrix.shape != (n, n) or not np.all(np.isfinite(matrix)):
            raise ValueError(f"{name} must be a finite {n} x {n} covariance")
        if not np.allclose(matrix, matrix.T):
            raise ValueError(f"{name} must be symmetric")
        np.linalg.cholesky(matrix)
    if cfg.max_iteration < 1 or cfg.window_size < 1:
        raise ValueError("Iteration count and window size must be positive")
    if cfg.robust_delta <= 0 or cfg.threshold_iteration < 0:
        raise ValueError("Invalid robust delta or convergence threshold")


def workspace_root():
    """The self-contained workspace (or installed package), never the parent repo."""
    return Path(__file__).resolve().parents[1]


def resource_path(name):
    return workspace_root() / "config" / name


def resolve_data_path(config):
    path = Path(config.data_path).expanduser()
    if path.is_absolute():
        return path
    # Relative to the distribution, never to an unrelated process cwd.
    return workspace_root() / "data" / path
