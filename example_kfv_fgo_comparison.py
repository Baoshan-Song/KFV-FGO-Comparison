from experiment_common import make_config, make_data, run_and_display


def main():
    data = make_data(anchor_radius=105)
    # Reproduce the MATLAB recurrence; see the mathematical caveat in README.
    config = make_config(kfv_mode="riEKF", max_iteration=10,
                         robust_kernel="huber", window_size=1, imitate_kfv=True)
    first = {"name": "KFV", "kind": "kfv", "config": config}
    second = {"name": "Re-FGO", "kind": "fgo", "config": config}
    return run_and_display(data, first, second)


if __name__ == "__main__":
    main()
