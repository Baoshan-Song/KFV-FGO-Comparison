from pathlib import Path
from config.config import load_config


def load_settings(path=None):
    return load_config(path or Path(__file__).with_name("swfgo_test.json"))
