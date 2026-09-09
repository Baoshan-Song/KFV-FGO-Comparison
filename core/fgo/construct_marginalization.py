def construct_marginalization(estimator, states_to_remove):
    return estimator.marginalize(states_to_remove)
