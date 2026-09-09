def estimate_fgo(estimator):
	"""Run Gauss-Newton on a core.fgo.factor_graph.FactorGraph."""
	return estimator.estimate()

__all__ = ["estimate_fgo"]
