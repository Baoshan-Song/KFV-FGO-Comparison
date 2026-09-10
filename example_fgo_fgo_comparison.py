from experiment_common import make_config, make_data, run_and_display


def main(): 
  data = make_data(anchor_radius=105)
  first = {"name": "FGO-standard", "kind": "fgo",
           "config": make_config(robust_kernel="none", imitate_kfv=True)}
  second = {"name": "FGO-robust", "kind": "fgo",
            "config": make_config(robust_kernel="none", robust_delta=2.0,
                                   imitate_kfv=False,window_size=20)}
  return run_and_display(data, first, second)


if __name__ == "__main__":
  main()