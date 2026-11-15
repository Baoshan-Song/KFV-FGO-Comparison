import json
from pathlib import Path

def test_primary_notebook_compiles_and_redirects():
    root=Path(__file__).resolve().parents[2]
    notebook=json.loads((root/'KFV-FGO-Colab.ipynb').read_text())
    for cell in notebook['cells']:
        if cell['cell_type']=='code':compile(''.join(cell['source']),'colab','exec')
    # The three workspaces each expose both data sources. Article reproduction
    # is a separate script, not a fourth notebook module.
    code_cells=[cell for cell in notebook['cells'] if cell['cell_type']=='code']
    assert len(code_cells) == 4
    for cell, module in zip(code_cells[1:], ['KFV vs KFV','KFV vs FGO','FGO vs FGO']):
        assert f'display(comparison_ui("{module}"))' in ''.join(cell['source'])
    setup=''.join(code_cells[0]['source'])
    assert '"notebook_transport"' in setup and 'close_panels()' in setup
    redirect=json.loads((root/'kfv_fgo/Colab.ipynb').read_text())
    assert all(c['cell_type']=='markdown' for c in redirect['cells'])
