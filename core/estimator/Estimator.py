from abc import ABC, abstractmethod


class Estimator(ABC):
	def __init__(self, config, data):
		self.config = config
		self.data = data
		self.results = {}

	@abstractmethod
	def run(self):
		raise NotImplementedError

__all__ = ["Estimator"]
