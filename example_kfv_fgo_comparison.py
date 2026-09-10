from experiment_common import make_config, make_data, run_and_display


def main():
    data = make_data(anchor_radius=105)
    config = make_config(kfv_mode = "riEKF",max_iteration =20)
    first = {"name": "KFV", "kind": "kfv", "config": config}
    second = {"name": "Re-FGO", "kind": "fgo",
                        "config": make_config(imitate_kfv=True,kfv_mode = "riEKF",max_iteration =20)}
    return run_and_display(data, first, second)


if __name__ == "__main__":
    main()
