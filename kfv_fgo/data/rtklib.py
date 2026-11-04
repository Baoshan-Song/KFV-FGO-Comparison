"""pyRTKLIB adapter reproducing the MatRTKLIB Gobs/Gsat L1 pipeline.

All public units are metres, seconds, Hz, dB-Hz; angles in this module are radians.
Native RTKLIB objects are confined here. Estimators receive immutable epoch packets.
"""

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ..model.gnss_ins import ecef_to_llh


def rtklib():
    try:
        import pyrtklib
    except ImportError as exc:
        raise ImportError(
            'Real data requires pip install "kfv-fgo-comparison[real]"'
        ) from exc
    return pyrtklib


def native_array(r, values):
    values = np.asarray(values, dtype=float).reshape(-1)
    out = r.Arr1Ddouble(len(values))
    for i, v in enumerate(values):
        out[i] = float(v)
    return out


@dataclass(frozen=True)
class GnssEpoch:
    index: int
    gps_unix: float
    sat_ids: np.ndarray
    positions: np.ndarray
    clocks: np.ndarray
    pseudorange: np.ndarray
    snr: np.ndarray
    frequencies: np.ndarray
    tgd: np.ndarray


class RtklibBackend:
    def __init__(self, obs_file, nav_file):
        r = self.r = rtklib()
        obs, nav, sta = r.obs_t(), r.nav_t(), r.sta_t()
        self.epochs = []
        try:
            for file in (obs_file, nav_file):
                path = Path(file).resolve()
                if not path.is_file():
                    raise FileNotFoundError(path)
                # Native Windows RTKLIB requires absolute paths using backslashes.
                if r.readrnx(str(path), 1, "", obs, nav, sta) <= 0:
                    raise ValueError(f"Cannot read RINEX: {path}")
            r.uniqnav(nav)
            self.ion = np.array([nav.ion_gps[i] for i in range(8)])
            tgd = {}
            for i in range(nav.n):
                eph = nav.eph[i]
                tgd.setdefault(eph.sat, float(eph.tgd[0]) * r.CLIGHT)
            start = 0
            while start < obs.n:
                first = obs.data[start]
                end = start + 1
                while (
                    end < obs.n
                    and r.timediff(obs.data[end].time, first.time) <= r.DTTOL
                ):
                    end += 1
                count = end - start
                rs = r.Arr1Ddouble(6 * count)
                clocks = r.Arr1Ddouble(2 * count)
                variance = r.Arr1Ddouble(count)
                health = r.Arr1Dint(count)
                r.satposs(
                    first.time,
                    first,
                    count,
                    nav,
                    r.EPHOPT_BRDC,
                    rs,
                    clocks,
                    variance,
                    health,
                )
                rows = []
                for j in range(count):
                    item = obs.data[start + j]
                    prn = r.Arr1Dint(1)
                    system = r.satsys(item.sat, prn)
                    position = np.array([rs[j * 6 + k] for k in range(3)])
                    p = float(item.P[0])
                    frequency = r.sat2freq(item.sat, item.code[0], nav)
                    if p == 0 or frequency <= 0 or np.linalg.norm(position) == 0:
                        continue
                    # Match Gsat: retain QZSS regardless of health flag.
                    if health[j] != 0 and system != r.SYS_QZS:
                        continue
                    snr = float(item.SNR[0]) * r.SNR_UNIT
                    if snr == 0:
                        snr = np.nan
                    rows.append(
                        (
                            item.sat,
                            position,
                            float(clocks[j * 2]) * r.CLIGHT,
                            p,
                            snr,
                            frequency,
                            tgd.get(item.sat, 0.0),
                        )
                    )
                rows.sort(key=lambda x: x[0])
                if not rows:
                    raise ValueError(
                        f"GNSS epoch {len(self.epochs) + 1}: no valid satellites"
                    )
                arrays = [np.array([row[k] for row in rows]) for k in range(7)]
                arrays[0] = arrays[0].astype(int)
                for array in arrays:
                    array.flags.writeable = False
                self.epochs.append(
                    GnssEpoch(
                        len(self.epochs), first.time.time + first.time.sec, *arrays
                    )
                )
                start = end
            if not self.epochs:
                raise ValueError("RINEX has no GNSS epochs")
        finally:
            r.freeobs(obs)
            r.freenav(nav, 0xFF)
        self.metadata = {
            "library": "pyrtklib",
            "version": "0.2.7",
            "rtklib": f"{r.VER_RTKLIB} {r.PATCH_LEVEL}",
        }

    def geometry(self, position, epoch):
        r = self.r
        native_position = native_array(r, position)
        llh = native_array(r, ecef_to_llh(position))
        ion = native_array(r, self.ion)
        time = r.gtime_t()
        time.time = int(epoch.gps_unix)
        time.sec = float(epoch.gps_unix - int(epoch.gps_unix))
        distances = []
        lines = []
        elevations = []
        ions = []
        trops = []
        for satellite, frequency in zip(epoch.positions, epoch.frequencies):
            los = r.Arr1Ddouble(3)
            azel = r.Arr1Ddouble(2)
            distance = r.geodist(native_array(r, satellite), native_position, los)
            r.satazel(llh, los, azel)
            distances.append(distance)
            lines.append([los[i] for i in range(3)])
            elevations.append(azel[1])
            ions.append(r.ionmodel(time, ion, llh, azel) * (r.FREQ1 / frequency) ** 2)
            trops.append(r.tropmodel(time, llh, azel, 0.7))
        return {
            "range": np.array(distances),
            "los": np.array(lines),
            "elevation": np.array(elevations),
            "ionosphere": np.array(ions),
            "troposphere": np.array(trops),
        }
