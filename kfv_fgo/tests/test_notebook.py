import io
import json
import sys
import types
import zipfile
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def code_cells():
    notebook = json.loads((ROOT / "Colab.ipynb").read_text(encoding="utf-8"))
    return [
        "".join(cell["source"])
        for cell in notebook["cells"]
        if cell["cell_type"] == "code"
    ]


def test_notebook_code_compiles_and_existing_workspace_is_preserved():
    for source in code_cells():
        compile(source, "Colab.ipynb", "exec")
    setup = code_cells()[0].replace(
        'SOURCE_MODE = "upload"', 'SOURCE_MODE = "existing"'
    )
    setup = setup.replace('WORKSPACE_PATH = ""', f"WORKSPACE_PATH = {str(ROOT)!r}")
    before = (ROOT / "core/estimator.py").read_bytes()
    ns = {}
    exec(setup, ns)
    exec(setup, ns)
    assert ns["WORKSPACE"] == ROOT
    assert (ROOT / "core/estimator.py").read_bytes() == before


@pytest.mark.parametrize(
    "path", ["../escape.py", "/absolute.py", "C:/escape.py", "..\\escape.py"]
)
def test_notebook_upload_rejects_escaping_paths(monkeypatch, tmp_path, path):
    archive = io.BytesIO()
    with zipfile.ZipFile(archive, "w") as z:
        z.writestr(path, "bad")
    google = types.ModuleType("google")
    module = types.ModuleType("google.colab")
    module.files = types.SimpleNamespace(
        upload=lambda: {"bundle.zip": archive.getvalue()}
    )
    google.colab = module
    monkeypatch.setitem(sys.modules, "google", google)
    monkeypatch.setitem(sys.modules, "google.colab", module)
    import tempfile

    monkeypatch.setattr(tempfile, "mkdtemp", lambda **kwargs: str(tmp_path))
    with pytest.raises(ValueError, match="unsafe path"):
        exec(code_cells()[0], {})


@pytest.mark.parametrize("prefix", ["", "kfv_fgo/", "repository-python_colab/kfv_fgo/"])
def test_notebook_accepts_source_zip_layouts(monkeypatch, tmp_path, prefix):
    archive = io.BytesIO()
    with zipfile.ZipFile(archive, "w") as z:
        for filename in (
            "pyproject.toml",
            "core/estimator.py",
            "examples/compare_matlab_results.py",
        ):
            z.writestr(prefix + filename, "")
    google = types.ModuleType("google")
    module = types.ModuleType("google.colab")
    module.files = types.SimpleNamespace(
        upload=lambda: {"source.zip": archive.getvalue()}
    )
    google.colab = module
    monkeypatch.setitem(sys.modules, "google", google)
    monkeypatch.setitem(sys.modules, "google.colab", module)
    import tempfile

    monkeypatch.setattr(tempfile, "mkdtemp", lambda **kwargs: str(tmp_path))
    ns = {}
    exec(code_cells()[0], ns)
    assert ns["WORKSPACE"] == tmp_path / prefix
