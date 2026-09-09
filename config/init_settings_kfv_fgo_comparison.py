from pathlib import Path
from config import load_config


def load_settings(path=None):
    return load_config(path or Path(__file__).with_name("kfv_fgo_comparison.json"))
