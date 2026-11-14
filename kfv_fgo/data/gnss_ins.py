"""RINEX/Xsens data loading, MATLAB-compatible GNSS/IMU alignment."""

from time import perf_counter
from collections import OrderedDict
from threading import RLock

import numpy as np
from scipy.interpolate import interp1d

from ..config.settings import resolve_data_path
from ..model.gnss_ins import GnssObservationModel, ImuPropagationModel
from .rtklib import RtklibBackend


class ImuData:
    def __init__(self, path):
        if str(path).endswith('.npz'):
            with np.load(path, allow_pickle=False) as cached:
                values = np.column_stack((cached['time_s'],cached['quat_xyzw'],cached['acc_mps2']))
            self.time = values[:,0]
        else:
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
        # Bracket in O(log N), rather than scanning the entire IMU recording
        # for both ends of every GNSS interval.
        second = min(int(np.searchsorted(self.time, t, side='right')), len(self.time)-1)
        first = max(second-1, 0) if self.time[second] > t else second
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
        indices = slice(int(np.searchsorted(self.time,start,side='left')),
                        int(np.searchsorted(self.time,end,side='right')))
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
    _loaded = OrderedDict()
    _lock = RLock()

    @classmethod
    def load(cls, config, progress=None, stop=None):
        from run_cache import input_fingerprint
        import json
        key=(input_fingerprint(config), config.leap_seconds, json.dumps(config.gnss,sort_keys=True),
             config.imu.get('gravity_norm',9.81))
        with cls._lock:
            if key in cls._loaded:
                if progress: progress('Reusing parsed GNSS/IMU observations')
                cls._loaded.move_to_end(key)
                return cls._loaded[key]
        data=cls(config,progress=progress,stop=stop)
        with cls._lock:
            cls._loaded[key]=data
            while len(cls._loaded)>2:cls._loaded.popitem(last=False)
        return data

    def __init__(self, config, progress=None, stop=None):
        started = perf_counter()
        def report(message):
            if stop is not None and stop.is_set():raise InterruptedError('Stopped while loading observations')
            if progress:progress(message)
        report('Checking GNSS/IMU files')
        self.path = resolve_data_path(config)
        for name in ("f9p_navi.obs", "brdm.rnx"):
            if not (self.path / name).is_file():
                raise FileNotFoundError(self.path / name)
        imu_path = self.path / 'xsens_imu.csv'
        if not imu_path.is_file():
            imu_path = self.path / 'xsens_imu.npz'
        if not imu_path.is_file():
            raise FileNotFoundError(f'{self.path}: need xsens_imu.csv or xsens_imu.npz')
        self.backend = RtklibBackend(self.path / "f9p_navi.obs", self.path / "brdm.rnx", progress=report)
        self.epochs = self.backend.epochs
        self.timestamps = np.array(
            [ep.gps_unix - config.leap_seconds for ep in self.epochs]
        )
        self.num_steps = len(self.timestamps)
        report('Reading IMU samples')
        imu = ImuData(imu_path)
        self.batches = []
        for k,(a,b) in enumerate(zip(self.timestamps[:-1],self.timestamps[1:])):
            if k%100==0:report(f'Aligning IMU intervals {k}/{self.num_steps-1}')
            self.batches.append(imu.interval(a,b))
        report('Computing initial GNSS state')
        self.observation_model = GnssObservationModel(self.backend, config.gnss)
        self.initial_state = self.observation_model.initialize(self.epochs[0])
        self.gravity_norm = config.imu.get("gravity_norm", 9.81)
        self.load_seconds = perf_counter() - started
        report(f'Ready: {self.num_steps-1} GNSS updates ({self.load_seconds:.1f} s loading)')

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
