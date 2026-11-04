"""RINEX/Xsens data loading, MATLAB-compatible GNSS/IMU alignment."""

from time import perf_counter

import numpy as np
from scipy.interpolate import interp1d

from ..config.settings import resolve_data_path
from ..model.gnss_ins import GnssObservationModel, ImuPropagationModel
from .rtklib import RtklibBackend


class ImuData:
    def __init__(self, path):
        values = np.loadtxt(
            path, delimiter=",", skiprows=1, usecols=(0, 4, 5, 6, 7, 29, 30, 31)
        )
        self.time = values[:, 0] * 1e-9
        _, indices = np.unique(self.time, return_index=True)
        indices.sort()
        self.time = self.time[indices]
        self.quat = values[indices, 1:5]
        self.acc = values[indices, 5:8]
        if (
            len(self.time) < 2
            or not np.all(np.isfinite(self.time))
            or np.any(np.diff(self.time) <= 0)
        ):
            raise ValueError(
                "IMU timestamps must be finite and increasing after duplicate removal"
            )
        self.interpolate_acc = interp1d(
            self.time, self.acc, axis=0, bounds_error=False, fill_value=np.nan
        )

    def quaternion_at(self, t):
        i = int(np.argmin(np.abs(self.time - t)))
        first, second = (
            (i, min(i + 1, len(self.time) - 1))
            if self.time[i] <= t
            else (max(i - 1, 0), i)
        )
        if first == second:
            return self.quat[first].copy()
        fraction = (t - self.time[first]) / (self.time[second] - self.time[first])
        q1, q2 = self.quat[first].copy(), self.quat[second].copy()
        dot = np.dot(q1, q2)
        if dot < 0:
            q2 = -q2
            dot = -dot
        if dot > 0.9995:
            q = (1 - fraction) * q1 + fraction * q2
            return q / np.linalg.norm(q)
        theta = np.arccos(np.clip(dot, -1, 1))
        scaled = theta * fraction
        return (np.cos(scaled) - dot * np.sin(scaled) / np.sin(theta)) * q1 + np.sin(
            scaled
        ) / np.sin(theta) * q2

    def interval(self, start, end):
        if end <= start:
            raise ValueError("GNSS epochs must be strictly increasing")
        indices = np.flatnonzero((self.time >= start) & (self.time <= end))
        return {
            "time": np.r_[start, self.time[indices], end],
            "acc": np.vstack(
                (
                    self.interpolate_acc(start),
                    self.acc[indices],
                    self.interpolate_acc(end),
                )
            ),
            "quat": np.vstack(
                (self.quaternion_at(start), self.quat[indices], self.quaternion_at(end))
            ),
        }


class GnssImuDataset:
    def __init__(self, config):
        started = perf_counter()
        self.path = resolve_data_path(config)
        for name in ("f9p_navi.obs", "brdm.rnx", "xsens_imu.csv"):
            if not (self.path / name).is_file():
                raise FileNotFoundError(self.path / name)
        self.backend = RtklibBackend(self.path / "f9p_navi.obs", self.path / "brdm.rnx")
        self.epochs = self.backend.epochs
        self.timestamps = np.array(
            [ep.gps_unix - config.leap_seconds for ep in self.epochs]
        )
        self.num_steps = len(self.timestamps)
        imu = ImuData(self.path / "xsens_imu.csv")
        self.batches = [
            imu.interval(a, b)
            for a, b in zip(self.timestamps[:-1], self.timestamps[1:])
        ]
        self.observation_model = GnssObservationModel(self.backend, config.gnss)
        self.initial_state = self.observation_model.initialize(self.epochs[0])
        self.gravity_norm = config.imu.get("gravity_norm", 9.81)
        self.load_seconds = perf_counter() - started

    def get_propagation_input(self, index):
        return self.batches[index - 1]

    def get_measurement(self, index):
        return self.epochs[index]

    def propagate(self, x, index, cfg):
        return ImuPropagationModel(cfg.q, self.gravity_norm).propagate(
            x, self.get_propagation_input(index)
        )

    def linearize(self, x, index):
        return self.observation_model.linearize(x, self.get_measurement(index))
