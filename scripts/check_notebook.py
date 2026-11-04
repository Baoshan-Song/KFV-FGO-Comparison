"""Execute every comparison cell locally, including the real-data GUI."""
import json
from pathlib import Path
from unittest.mock import patch
from IPython.display import HTML
import nbformat
nbformat.validate(nbformat.read('KFV-FGO-Colab.ipynb',as_version=4))
nb=json.loads(Path('KFV-FGO-Colab.ipynb').read_text())
ns={};previews=[]
for i,cell in enumerate(nb['cells'][1:],1):
    if cell['cell_type']!='code':continue
    with patch('IPython.display.display',side_effect=lambda v: previews.append(v.data) if isinstance(v,HTML) else None):
        exec(compile(''.join(cell['source']),f'colab_cell_{i}','exec'),ns)
    print('PASS notebook cell',i)
Path('artifacts').mkdir(exist_ok=True)
for name,content in [('colab-preview',previews[1]),('real-colab-preview',previews[-1])]:
    Path(f'artifacts/{name}.html').write_text('<!doctype html><meta charset="utf-8">'+content)
print('HTML previews written')
