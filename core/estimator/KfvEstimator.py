import numpy as np
from .Estimator import Estimator
from ..filter import ekf, miekf, rekf, rmiekf
from config.config import convert_kfv_to_fgo, motion, motion_jacobian, range_measurement, range_jacobian


class KfvEstimator(Estimator):
	def run(self):
		cfg = self.config.kfv
		count = self.data["num_steps"]
		emitters = self.data["emitter_positions"]
		state = np.r_[self.data["true_positions"][:, 0], self.data["true_velocities"][:, 0]] + cfg.err_x0
		covariance = cfg.p0.copy()
		estimates = np.zeros((4, count))
		estimates[:, 0] = state
		debug_info = [None] * count
		methods = {"EKF": ekf, "iEKF": miekf, "rEKF": rekf, "riEKF": rmiekf}
		if cfg.mode not in methods:
			raise ValueError(f"Unsupported KFV mode: {cfg.mode}")
		for index in range(1, count):
			arguments = (state, covariance, cfg.dt, cfg.omega, motion, motion_jacobian, cfg.q,
						 self.data["toa_measurements"][:, index], emitters, range_measurement,
						 range_jacobian, cfg.r)
			if cfg.mode == "EKF": result = methods[cfg.mode](*arguments)
			elif cfg.mode == "iEKF": result = methods[cfg.mode](*arguments, cfg.max_iteration, cfg.threshold_iteration)
			elif cfg.mode == "rEKF": result = methods[cfg.mode](*arguments, cfg.robust_kernel, cfg.robust_delta)
			else: result = methods[cfg.mode](*arguments, cfg.max_iteration, cfg.threshold_iteration,
											 cfg.robust_kernel, cfg.robust_delta)
			state, covariance, _, _, debug_info[index] = result
			estimates[:, index] = state
		return {"X": estimates, "debug_info": debug_info}

	def convert_kfv_config_to_fgo(self):
		return convert_kfv_to_fgo(self.config)

__all__ = ["KfvEstimator"]
