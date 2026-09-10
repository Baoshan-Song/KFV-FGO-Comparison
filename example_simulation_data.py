import argparse
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from data.circle_eval import generate_data
from config.config import range_jacobian


def compute_nonlinearity_metrics(data):
  """Computes Frobenius norm of J^T J along the trajectory to quantify geometric non-linearity."""
  true_positions = data["true_positions"]
  emitter_positions = data["emitter_positions"]
  num_steps = true_positions.shape[1]

  hessian_norms = np.zeros(num_steps)
  for k in range(num_steps):
    x_k = np.array(
        [true_positions[0, k], true_positions[1, k], 0.0, 0.0]
    )  # [x, y, vx, vy]
    jacobians = np.vstack(
        [range_jacobian(x_k, emitter) for emitter in emitter_positions.T]
    )
    hessian_norms[k] = np.linalg.norm(jacobians.T @ jacobians, "fro")

  return hessian_norms


def plot_diagnostics(data, anchor_radius, gmm_params, save_fig_path=None):
  """Plots 2 figures evaluating Non-linearity (Geometry/Hessian) and Non-Gaussianity (Residual GMM distribution)."""
  true_pos = data["true_positions"]
  emitters = data["emitter_positions"]
  toa_meas = data["toa_measurements"]
  num_steps = data["num_steps"]

  # Compute true range residuals (Measurement minus True Geometry Range)
  true_ranges = np.zeros_like(toa_meas)
  for i in range(emitters.shape[1]):
    diff = true_pos - emitters[:, i : i + 1]
    true_ranges[i, :] = np.linalg.norm(diff, axis=0)
  residuals = (toa_meas - true_ranges).flatten()

  hessian_norms = compute_nonlinearity_metrics(data)

  fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5), dpi=100)

  # --------------------------------------------------------------------------
  # Subplot 1: Geometric Non-Linearity (Trajectory & Hessian Norm Heatmap)
  # --------------------------------------------------------------------------
  sc = ax1.scatter(
      true_pos[0, :],
      true_pos[1, :],
      c=hessian_norms,
      cmap="plasma",
      s=15,
      zorder=3,
      label="UAV Path (Color = $\|\mathbf{J}^T \mathbf{J}\|_F$)",
  )
  cbar = plt.colorbar(sc, ax=ax1)
  cbar.set_label("Hessian Matrix Norm $\|\mathbf{J}^T \mathbf{J}\|_F$ (Non-linearity)")

  # Plot UWB Emitters (Anchors)
  ax1.scatter(
      emitters[0, :],
      emitters[1, :],
      c="red",
      marker="^",
      s=120,
      edgecolors="black",
      zorder=4,
      label=f"UWB Anchors (R = {anchor_radius}m)",
  )

  # Draw Anchor Array Circle Reference
  circle = plt.Circle(
      (0, 0),
      anchor_radius,
      color="red",
      fill=False,
      linestyle="--",
      alpha=0.5,
      label="Anchor Array Boundary",
  )
  ax1.add_patch(circle)

  mean_norm = np.mean(hessian_norms)
  ax1.set_title(
        f"1. Geometric Non-Linearity\nAnchor Radius R = {anchor_radius}m "
        f"(Mean $\|\\mathbf{{J}}^T \\mathbf{{J}}\\|_F$ = {mean_norm:.2f})"
    )
  ax1.set_xlabel("X Position [m]")
  ax1.set_ylabel("Y Position [m]")
  ax1.grid(True, linestyle=":", alpha=0.6)
  ax1.legend(loc="upper right", fontsize=8)
  ax1.axis("equal")

  # --------------------------------------------------------------------------
  # Subplot 2: Measurement Non-Gaussianity (GMM Residual PDF vs Normal)
  # --------------------------------------------------------------------------
  w1, w2, mu1, mu2, sig1, sig2 = gmm_params

  # Plot Empirical Residual Histogram
  count, bins, _ = ax2.hist(
      residuals,
      bins=60,
      density=True,
      alpha=0.5,
      color="skyblue",
      edgecolor="navy",
      label="Simulated Residuals (GMM)",
  )

  # Plot Analytical GMM PDF Curve
  x_pdf = np.linspace(np.min(residuals) - 2, np.max(residuals) + 2, 500)
  pdf1 = (1 / (sig1 * np.sqrt(2 * np.pi))) * np.exp(
      -0.5 * ((x_pdf - mu1) / sig1) ** 2
  )
  pdf2 = (1 / (sig2 * np.sqrt(2 * np.pi))) * np.exp(
      -0.5 * ((x_pdf - mu2) / sig2) ** 2
  )
  gmm_pdf = w1 * pdf1 + w2 * pdf2
  ax2.plot(
      x_pdf,
      gmm_pdf,
      "r-",
      lw=2.5,
      label=(
          f"True GMM PDF (w2={w2:.2f}, μ2={mu2}m)\nComponent 1: N({mu1},"
          f" {sig1}²)\nComponent 2: N({mu2}, {sig2}²)"
      ),
  )

  # Overlay Gaussian Assumption Curve (How EKF views the noise)
  nominal_pdf = (1 / (sig1 * np.sqrt(2 * np.pi))) * np.exp(
      -0.5 * (x_pdf / sig1) ** 2
  )
  ax2.plot(
      x_pdf,
      nominal_pdf,
      "k--",
      lw=1.5,
      label=f"Standard Gaussian Assumption N(0, {sig1}²)",
  )

  ax2.set_title(
      f"2. Measurement Non-Gaussianity (NLOS Heavy-Tail)\nMean Offset ="
      f" {np.mean(residuals):.2f}m, Std = {np.std(residuals):.2f}m"
  )
  ax2.set_xlabel("UWB Range Measurement Residual [m]")
  ax2.set_ylabel("Probability Density")
  ax2.grid(True, linestyle=":", alpha=0.6)
  ax2.legend(loc="upper right", fontsize=8)

  plt.tight_layout()

  if save_fig_path:
    plt.savefig(save_fig_path, dpi=300)
    print(f"📊 Diagnostic plots saved to: {save_fig_path}")

  plt.show()


def generate_custom_sim_data(
    anchor_radius: float = 200.0,
    gmm_w1: float = 0.8,
    gmm_w2: float = 0.2,
    gmm_mu1: float = 0.0,
    gmm_mu2: float = 30.0,
    gmm_sigma1: float = 0.1,
    gmm_sigma2: float = 5.0,
    save_filename: str = "circle_cv_gmm_L4.mat",
    plot_result: bool = True,
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

  print("✅ Simulation Data Generated & Saved Successfully!")
  print(f"   • File Path        : {file_path}")
  print(f"   • Anchor Radius    : {anchor_radius} m")
  print(f"   • GMM Weights (w)  : ({gmm_w1:.2f}, {gmm_w2:.2f})")
  print(f"   • GMM Means (μ)    : ({gmm_mu1}, {gmm_mu2}) m")
  print(f"   • GMM Sigmas (σ)   : ({gmm_sigma1}, {gmm_sigma2}) m")

  if plot_result:
    gmm_params = (gmm_w1, gmm_w2, gmm_mu1, gmm_mu2, gmm_sigma1, gmm_sigma2)
    fig_save_path = data_dir / (file_path.stem + "_diagnostics.png")
    plot_diagnostics(
        data, anchor_radius, gmm_params, save_fig_path=fig_save_path
    )

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
  parser.add_argument(
      "--gmm_w1", type=float, default=0.8, help="Weight of Component 1"
  )
  parser.add_argument(
      "--gmm_w2", type=float, default=0.2, help="Weight of Component 2"
  )
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
  parser.add_argument(
      "--plot",
      action="store_true",
      default=True,
      help="Generate diagnostic plots",
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
      plot_result=args.plot,
  )
  
  
# import argparse
# from pathlib import Path
# import numpy as np
# from data.circle_eval import generate_data


# def generate_custom_sim_data(
#     anchor_radius: float = 200.0,
#     gmm_w1: float = 0.8,
#     gmm_w2: float = 0.2,
#     gmm_mu1: float = 0.0,
#     gmm_mu2: float = 30.0,
#     gmm_sigma1: float = 0.1,
#     gmm_sigma2: float = 5.0,
#     save_filename: str = "circle_cv_gmm_L4.mat",
# ):
#     root = (
#         Path(__file__).parent.parent
#         if Path(__file__).parent.name == "data"
#         else Path(__file__).parent
#     )
#     data_dir = root / "data"
#     data_dir.mkdir(parents=True, exist_ok=True)
#     file_path = data_dir / save_filename

#     gmm_weights = (gmm_w1, gmm_w2)
#     gmm_means = (gmm_mu1, gmm_mu2)
#     gmm_sigmas = (gmm_sigma1, gmm_sigma2)

#     data = generate_data(
#         emitter_radius=anchor_radius,
#         gmm_weights=gmm_weights,
#         gmm_means=gmm_means,
#         gmm_sigmas=gmm_sigmas,
#         save_path=file_path,
#     )

#     print(f"✅ Simulation Data Generated & Saved Successfully!")
#     print(f"   • File Path        : {file_path}")
#     print(f"   • Anchor Radius    : {anchor_radius} m")
#     print(f"   • GMM Weights (w)  : ({gmm_w1:.2f}, {gmm_w2:.2f})")
#     print(f"   • GMM Means (μ)    : ({gmm_mu1}, {gmm_mu2}) m")
#     print(f"   • GMM Sigmas (σ)   : ({gmm_sigma1}, {gmm_sigma2}) m")

#     return str(file_path)


# if __name__ == "__main__":
#     parser = argparse.ArgumentParser(
#         description="Generate synthetic simulation data."
#     )
#     parser.add_argument(
#         "--anchor_radius",
#         type=float,
#         default=200.0,
#         help="Anchor radius in meters",
#     )
#     parser.add_argument("--gmm_w1", type=float, default=0.8, help="Weight of Component 1")
#     parser.add_argument("--gmm_w2", type=float, default=0.2, help="Weight of Component 2")
#     parser.add_argument(
#         "--gmm_mu1", type=float, default=0.0, help="Mean of Component 1"
#     )
#     parser.add_argument(
#         "--gmm_mu2", type=float, default=30.0, help="Mean of Component 2"
#     )
#     parser.add_argument(
#         "--gmm_sigma1",
#         type=float,
#         default=0.1,
#         help="Standard deviation of Component 1",
#     )
#     parser.add_argument(
#         "--gmm_sigma2",
#         type=float,
#         default=5.0,
#         help="Standard deviation of Component 2",
#     )
#     parser.add_argument(
#         "--output_filename",
#         type=str,
#         default="circle_cv_gmm_L4.mat",
#         help="Output mat filename",
#     )

#     args = parser.parse_args()

#     generate_custom_sim_data(
#         anchor_radius=args.anchor_radius,
#         gmm_w1=args.gmm_w1,
#         gmm_w2=args.gmm_w2,
#         gmm_mu1=args.gmm_mu1,
#         gmm_mu2=args.gmm_mu2,
#         gmm_sigma1=args.gmm_sigma1,
#         gmm_sigma2=args.gmm_sigma2,
#         save_filename=args.output_filename,
#     )
# # from pathlib import Path
# # import numpy as np
# # from config.config import motion, motion_jacobian, range_measurement, range_jacobian
# # from data.circle_eval import generate_data, load_data
# # from kfv_fgo.filters import ekf


# # def main():
# #     root = Path(__file__).parent
# #     path = root / "data" / "circle_cv_gmm_L2.mat"
# #     data = load_data(path) if path.exists() else generate_data(gmm_weights=(1.0, 0.0), gmm_sigmas=(0.1, 10.0))
# #     dt = float(np.asarray(data["dt"]).squeeze())
# #     omega = float(np.asarray(data["omega"]).squeeze())
# #     x = np.array([200.0, 0.0, 0.0, omega]); covariance = np.diag([50.0, 50.0, 0.1, 0.1])
# #     estimates = np.zeros((4, data["num_steps"])); estimates[:, 0] = x
# #     norms = np.zeros(data["num_steps"])
# #     for index in range(1, data["num_steps"]):
# #         x, covariance, *_ = ekf(x, covariance, dt, omega, motion, motion_jacobian,
# #                                  np.diag([1e-2, 1e-2, 1e-4, 1e-4]), data["toa_measurements"][:, index],
# #                                  data["emitter_positions"], range_measurement, range_jacobian, 1.0)
# #         estimates[:, index] = x
# #         jacobians = np.vstack([range_jacobian(x, emitter) for emitter in data["emitter_positions"].T])
# #         norms[index] = np.linalg.norm(jacobians.T @ jacobians, "fro")
# #     error = np.linalg.norm(estimates[:2] - data["true_positions"], axis=0)
# #     print({"position_rmse": float(np.sqrt(np.mean(error ** 2))), "hessian_norms": norms})


# # if __name__ == "__main__": main()
