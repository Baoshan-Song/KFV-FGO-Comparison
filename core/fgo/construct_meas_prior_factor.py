def construct_meas_prior_factor(estimator, related_state):
    return estimator.mar_measurements(related_state)
