from __future__ import annotations

from pathlib import Path
import numpy as np
from scipy.io import loadmat, savemat


def _scalar(value):
    return float(np.asarray(value).squeeze())


def load_data(path: str | Path) -> dict:
    raw = loadmat(path)
    data = {key: value for key, value in raw.items() if not key.startswith("__")}
    for key in ("true_positions", "true_velocities", "emitter_positions", "toa_measurements"):
        if key in data:
            data[key] = np.asarray(data[key], dtype=float)
    data["num_steps"] = int(_scalar(data.get("num_steps", data["true_positions"].shape[1])))
    return data


def generate_data(num_steps=100, dt=1.0, radius=100.0, num_emitters=6,
                  emitter_radius=105.0, gmm_weights=(1.0, 0.0),
                  gmm_sigmas=(0.1, 10.0), seed=42, gmm_means=None,
                  save_path=None) -> dict:
    weights = np.asarray(gmm_weights, dtype=float)
    sigmas = np.asarray(gmm_sigmas, dtype=float)
    means = np.zeros_like(sigmas) if gmm_means is None else np.asarray(gmm_means, dtype=float)
    if weights.ndim != 1 or sigmas.ndim != 1 or means.ndim != 1:
        raise ValueError("GMM weights, means, and sigmas must be one-dimensional")
    if not (len(weights) == len(sigmas) == len(means)):
        raise ValueError("GMM weights, means, and sigmas must have the same length")
    if np.any(weights < 0) or not np.isfinite(weights).all() or weights.sum() <= 0:
        raise ValueError("GMM weights must be finite, non-negative, and have a positive sum")
    if np.any(sigmas < 0) or not np.isfinite(sigmas).all():
        raise ValueError("GMM sigmas must be finite and non-negative")

    rng = np.random.default_rng(seed)
    omega = 2 * np.pi / (num_steps * dt)
    angles = np.linspace(0, 2 * np.pi, num_steps)
    positions = radius * np.vstack((np.cos(angles), np.sin(angles)))
    velocities = omega * radius * np.vstack((-np.sin(angles), np.cos(angles)))
    emitter_angles = np.linspace(0, 2 * np.pi, num_emitters, endpoint=False)
    emitters = emitter_radius * np.vstack((np.cos(emitter_angles), np.sin(emitter_angles)))
    distances = np.linalg.norm(positions[:, None, :] - emitters[:, :, None], axis=0)
    components = rng.choice(len(weights), size=(num_emitters, num_steps), p=weights / weights.sum())
    noise = rng.normal(size=(num_emitters, num_steps)) * sigmas[components] + means[components]
    data = {"num_steps": num_steps, "dt": dt, "omega": omega,
            "true_positions": positions, "true_velocities": velocities,
            "emitter_positions": emitters, "toa_measurements": distances + noise}
    if save_path is not None:
        savemat(Path(save_path), data)
    return data
