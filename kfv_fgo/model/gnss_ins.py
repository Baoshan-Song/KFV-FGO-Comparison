"""Ten-state GNSS/INS propagation, pseudorange model and coordinates."""

import numpy as np


def ecef_to_llh(xyz):
    """RTKLIB's WGS84 conversion; latitude/longitude in radians."""
    r = np.asarray(xyz, dtype=float)
    e2 = (1 / 298.257223563) * (2 - 1 / 298.257223563)
    r2 = r[0] ** 2 + r[1] ** 2
    z, previous, v = r[2], 0.0, 6378137.0
    while abs(z - previous) >= 1e-4:
        previous = z
        sinp = z / np.sqrt(r2 + z * z)
        v = 6378137.0 / np.sqrt(1 - e2 * sinp * sinp)
        z = r[2] + v * e2 * sinp
    lat = (
        np.arctan(z / np.sqrt(r2))
        if r2 > 1e-12
        else (np.pi / 2 if r[2] > 0 else -np.pi / 2)
    )
    lon = np.arctan2(r[1], r[0]) if r2 > 1e-12 else 0.0
    return np.array([lat, lon, np.sqrt(r2 + z * z) - v])


def local_to_ecef(position):
    lat, lon, _ = ecef_to_llh(position)
    sl, cl, sp, cp = np.sin(lon), np.cos(lon), np.sin(lat), np.cos(lat)
    return np.array([[-sl, -sp * cl, cp * cl], [cl, -sp * sl, cp * sl], [0.0, cp, sp]])


def quaternion_rotation(q):
    q = np.asarray(q, dtype=float)
    norm = np.linalg.norm(q)
    if not np.isfinite(norm) or norm == 0:
        raise ValueError("Invalid zero/non-finite IMU quaternion")
    x, y, z, w = q / norm  # stored [x y z w]
    return np.array(
        [
            [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
            [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
            [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)],
        ]
    )


class ImuPropagationModel:
    def __init__(self, q, gravity_norm=9.81):
        self.q, self.gravity_norm = np.asarray(q), gravity_norm

    def propagate(self, x, batch):
        x = np.asarray(x)
        if x.shape != (10,):
            raise ValueError("IMU propagation needs a ten-state vector")
        dt = batch["time"][-1] - batch["time"][0]
        good = np.flatnonzero(
            np.all(np.isfinite(batch["acc"]), axis=1)
            & np.all(np.isfinite(batch["quat"]), axis=1)
        )
        if not len(good):
            raise ValueError("IMU interval has no valid acceleration/quaternion pair")
        first = good[0]
        ren = local_to_ecef(x[:3])
        reb = ren @ quaternion_rotation(batch["quat"][first])
        acceleration = reb @ (batch["acc"][first] - x[6:9]) + ren @ np.array(
            [0.0, 0.0, -self.gravity_norm]
        )
        predicted = x.copy()
        predicted[:3] += x[3:6] * dt
        predicted[3:6] += acceleration * dt
        F = np.eye(10)
        F[:3, 3:6] = np.eye(3) * dt
        F[3:6, 6:9] = -reb * dt
        return predicted, F, self.q


class GnssObservationModel:
    def __init__(self, backend, config):
        self.backend, self.config = backend, config

    def details(self, state, epoch):
        g = self.backend.geometry(state[:3], epoch)
        residual = (
            epoch.pseudorange
            - (g["range"] - epoch.clocks + g["ionosphere"] + g["troposphere"])
            - state[-1]
            - epoch.tgd
        )
        valid = np.isfinite(residual)
        if self.config.get("minimum_elevation_deg") is not None:
            valid &= g["elevation"] >= np.deg2rad(self.config["minimum_elevation_deg"])
        if self.config.get("minimum_snr") is not None:
            valid &= epoch.snr >= self.config["minimum_snr"]
        g.update(residual=residual, valid=valid)
        return g

    def linearize(self, state, epoch):
        g = self.details(state, epoch)
        valid = g["valid"]
        count = np.count_nonzero(valid)
        if count < 4:
            raise ValueError(
                f"GNSS epoch {epoch.index + 1}: need four valid satellites, found {count}"
            )
        H = np.zeros((count, len(state)))
        H[:, :3] = -g["los"][valid]
        H[:, -1] = 1
        elev = g["elevation"][valid]
        snr = epoch.snr[valid]
        tref = self.config.get("Tref", 45)
        a = self.config.get("a", 30)
        A = self.config.get("A", 30)
        fref = self.config.get("Fref", 10)
        snr = np.where(np.isfinite(snr), snr, tref)
        elev = np.where(np.isfinite(elev), elev, 1e-3)
        variance = (
            1
            / np.sin(np.maximum(elev, 1e-3)) ** 2
            * 10 ** (-(snr - tref) / a)
            * ((A * 10 ** ((fref - tref) / a) - 1) * (snr - tref) / (fref - tref) + 1)
        )
        if not np.all(np.isfinite(variance) & (variance > 0)):
            raise ValueError(
                f"GNSS epoch {epoch.index + 1}: invalid observation covariance"
            )
        return g["residual"][valid], H, np.diag(variance), epoch.sat_ids[valid]

    def initialize(self, epoch):
        state = np.zeros(4)
        for iteration in range(10):
            g = self.backend.geometry(state[:3], epoch)
            residual = (
                epoch.pseudorange
                - (g["range"] - epoch.clocks + g["ionosphere"] + g["troposphere"])
                - state[3]
                - epoch.tgd
            )
            good = np.isfinite(residual) & (g["elevation"] > np.deg2rad(10))
            if good.sum() < 4:
                raise ValueError(
                    f"SPP epoch {epoch.index + 1}: fewer than four satellites above 10 degrees"
                )
            H = np.column_stack((-g["los"][good], np.ones(good.sum())))
            root = np.sqrt(np.sin(g["elevation"][good]) / 0.5**2)
            delta = np.linalg.lstsq(
                root[:, None] * H, root * residual[good], rcond=None
            )[0]
            state += delta
            if np.linalg.norm(delta) < 1e-3:
                return np.r_[state[:3], np.zeros(6), state[3]]
        raise ValueError(
            f"SPP epoch {epoch.index + 1}: failed to converge within 10 iterations"
        )
