from experiment_common import make_config, make_data, run_and_display


def main():
    data = make_data(outlier_weight=0.2,anchor_radius=105)
    first = {"name": "EKF", "kind": "kfv",
                     "config": make_config(kfv_mode="EKF", robust_kernel="none")}
    second = {"name": "RIEKF", "kind": "kfv",
                        "config": make_config(kfv_mode="rEKF", max_iteration = 2)}
    return run_and_display(data, first, second)


if __name__ == "__main__":
    main()
