import argparse
from pathlib import Path
import numpy as np
from data.circle_eval import generate_data


def generate_custom_sim_data(
    anchor_radius: float = 200.0,
    gmm_w1: float = 0.8,
    gmm_w2: float = 0.2,
    gmm_mu1: float = 0.0,
    gmm_mu2: float = 30.0,
    gmm_sigma1: float = 0.1,
    gmm_sigma2: float = 5.0,
    save_filename: str = "circle_cv_gmm_L4.mat",
):
    root = (
        Path(__file__).parent.parent
        if Path(__file__).parent.name == "data"
        else Path(__file__).parent
    )
    data_dir = root / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    file_path = data_dir / save_filename

    gmm_weights = (gmm_w1, gmm_w2)
    gmm_means = (gmm_mu1, gmm_mu2)
    gmm_sigmas = (gmm_sigma1, gmm_sigma2)

    data = generate_data(
        emitter_radius=anchor_radius,
        gmm_weights=gmm_weights,
        gmm_means=gmm_means,
        gmm_sigmas=gmm_sigmas,
        save_path=file_path,
    )

    print(f"✅ Simulation Data Generated & Saved Successfully!")
    print(f"   • File Path        : {file_path}")
    print(f"   • Anchor Radius    : {anchor_radius} m")
    print(f"   • GMM Weights (w)  : ({gmm_w1:.2f}, {gmm_w2:.2f})")
    print(f"   • GMM Means (μ)    : ({gmm_mu1}, {gmm_mu2}) m")
    print(f"   • GMM Sigmas (σ)   : ({gmm_sigma1}, {gmm_sigma2}) m")

    return str(file_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate synthetic simulation data."
    )
    parser.add_argument(
        "--anchor_radius",
        type=float,
        default=200.0,
        help="Anchor radius in meters",
    )
    parser.add_argument("--gmm_w1", type=float, default=0.8, help="Weight of Component 1")
    parser.add_argument("--gmm_w2", type=float, default=0.2, help="Weight of Component 2")
    parser.add_argument(
        "--gmm_mu1", type=float, default=0.0, help="Mean of Component 1"
    )
    parser.add_argument(
        "--gmm_mu2", type=float, default=30.0, help="Mean of Component 2"
    )
    parser.add_argument(
        "--gmm_sigma1",
        type=float,
        default=0.1,
        help="Standard deviation of Component 1",
    )
    parser.add_argument(
        "--gmm_sigma2",
        type=float,
        default=5.0,
        help="Standard deviation of Component 2",
    )
    parser.add_argument(
        "--output_filename",
        type=str,
        default="circle_cv_gmm_L4.mat",
        help="Output mat filename",
    )

    args = parser.parse_args()

    generate_custom_sim_data(
        anchor_radius=args.anchor_radius,
        gmm_w1=args.gmm_w1,
        gmm_w2=args.gmm_w2,
        gmm_mu1=args.gmm_mu1,
        gmm_mu2=args.gmm_mu2,
        gmm_sigma1=args.gmm_sigma1,
        gmm_sigma2=args.gmm_sigma2,
        save_filename=args.output_filename,
    )
# from pathlib import Path
# import numpy as np
# from config.config import motion, motion_jacobian, range_measurement, range_jacobian
# from data.circle_eval import generate_data, load_data
# from kfv_fgo.filters import ekf


# def main():
#     root = Path(__file__).parent
#     path = root / "data" / "circle_cv_gmm_L2.mat"
#     data = load_data(path) if path.exists() else generate_data(gmm_weights=(1.0, 0.0), gmm_sigmas=(0.1, 10.0))
#     dt = float(np.asarray(data["dt"]).squeeze())
#     omega = float(np.asarray(data["omega"]).squeeze())
#     x = np.array([200.0, 0.0, 0.0, omega]); covariance = np.diag([50.0, 50.0, 0.1, 0.1])
#     estimates = np.zeros((4, data["num_steps"])); estimates[:, 0] = x
#     norms = np.zeros(data["num_steps"])
#     for index in range(1, data["num_steps"]):
#         x, covariance, *_ = ekf(x, covariance, dt, omega, motion, motion_jacobian,
#                                  np.diag([1e-2, 1e-2, 1e-4, 1e-4]), data["toa_measurements"][:, index],
#                                  data["emitter_positions"], range_measurement, range_jacobian, 1.0)
#         estimates[:, index] = x
#         jacobians = np.vstack([range_jacobian(x, emitter) for emitter in data["emitter_positions"].T])
#         norms[index] = np.linalg.norm(jacobians.T @ jacobians, "fro")
#     error = np.linalg.norm(estimates[:2] - data["true_positions"], axis=0)
#     print({"position_rmse": float(np.sqrt(np.mean(error ** 2))), "hessian_norms": norms})


# if __name__ == "__main__": main()
