import json
from pathlib import Path

def test_primary_notebook_compiles_and_redirects():
    root=Path(__file__).resolve().parents[2]
    notebook=json.loads((root/'KFV-FGO-Colab.ipynb').read_text())
    for cell in notebook['cells']:
        if cell['cell_type']=='code':compile(''.join(cell['source']),'colab','exec')
    last=''.join(notebook['cells'][-1]['source'])
    assert 'run_article_comparison' in last and 'save_real_experiment' in last
    pair=''.join(notebook['cells'][-2]['source'])
    assert 'run_real_pair' in pair and 'save_real_experiment' in pair
    redirect=json.loads((root/'kfv_fgo/Colab.ipynb').read_text())
    assert all(c['cell_type']=='markdown' for c in redirect['cells'])
