from abc import ABC, abstractmethod
from time import perf_counter
import numpy as np
from config.config import _config_to_dict


class Estimator(ABC):
	def __init__(self, config, data):
		self.config = config
		self.data = data
		self.results = {}

	def _simulation_output(self, result, started):
		result.update(solver=type(self).__name__, runtime_seconds=perf_counter()-started,
			timestamps=np.arange(self.data["num_steps"])*float(np.asarray(self.data.get("dt",1.)).squeeze()),
			config=_config_to_dict(self.config), backend={"model":"simulated ranges"})
		return result

	def run(self):
		stream = self.iter_epochs()
		while True:
			try: next(stream)
			except StopIteration as completed: return completed.value

	@abstractmethod
	def iter_epochs(self):
		raise NotImplementedError

__all__ = ["Estimator"]
