from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy.io import loadmat


def _scalar(value):
    return float(np.asarray(value).squeeze())


def load_data(path: str | Path) -> dict:
    raw = loadmat(path)
    data = {key: value for key, value in raw.items() if not key.startswith("__")}
    for key in (
        "true_positions",
        "true_velocities",
        "emitter_positions",
        "toa_measurements",
    ):
        if key in data:
            data[key] = np.asarray(data[key], dtype=float)
    data["num_steps"] = int(
        _scalar(data.get("num_steps", data["true_positions"].shape[1]))
    )
    return data


def generate_data(
    num_steps=100,
    dt=1.0,
    radius=100.0,
    num_emitters=6,
    emitter_radius=105.0,
    gmm_weights=(1.0, 0.0),
    gmm_sigmas=(0.1, 10.0),
    seed=42,
) -> dict:
    rng = np.random.default_rng(seed)
    omega = 2 * np.pi / (num_steps * dt)
    angles = np.linspace(0, 2 * np.pi, num_steps)
    positions = radius * np.vstack((np.cos(angles), np.sin(angles)))
    velocities = omega * radius * np.vstack((-np.sin(angles), np.cos(angles)))
    emitter_angles = np.linspace(0, 2 * np.pi, num_emitters, endpoint=False)
    emitters = emitter_radius * np.vstack(
        (np.cos(emitter_angles), np.sin(emitter_angles))
    )
    distances = np.linalg.norm(positions[:, None, :] - emitters[:, :, None], axis=0)
    components = rng.choice(
        len(gmm_weights),
        size=(num_emitters, num_steps),
        p=np.asarray(gmm_weights) / np.sum(gmm_weights),
    )
    noise = (
        rng.normal(size=(num_emitters, num_steps)) * np.asarray(gmm_sigmas)[components]
    )
    return {
        "num_steps": num_steps,
        "dt": dt,
        "omega": omega,
        "true_positions": positions,
        "true_velocities": velocities,
        "emitter_positions": emitters,
        "toa_measurements": distances + noise,
    }
