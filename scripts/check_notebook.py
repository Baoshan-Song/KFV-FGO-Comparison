import json
from pathlib import Path
from unittest.mock import patch
from IPython.display import HTML
import nbformat
nbformat.validate(nbformat.read("KFV-FGO-Colab.ipynb", as_version=4))

nb=json.loads(Path('KFV-FGO-Colab.ipynb').read_text())
ns={}
previews=[]
for i,cell in enumerate(nb['cells'][1:],1):
    with patch('IPython.display.display', side_effect=lambda v: previews.append(v.data) if isinstance(v,HTML) else None):
        exec(compile(''.join(cell['source']),f'colab_cell_{i}','exec'),ns)
    print('PASS notebook cell',i)
Path('artifacts').mkdir(exist_ok=True)
Path('artifacts/colab-preview.html').write_text('<!doctype html><meta charset="utf-8"><title>KFV / Re-FGO verification</title>'+previews[1])
print('HTML preview written')
