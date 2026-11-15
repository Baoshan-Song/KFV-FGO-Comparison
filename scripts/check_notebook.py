"""Check all three widget modules and smoke-test simulation/real routing."""
import argparse
from copy import copy
import json
from pathlib import Path
from unittest.mock import patch
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import nbformat
from IPython.display import HTML
from colab_comparison import prepare_comparison, MODULES
from experiment_common import run_pair, animate_pair
from real_experiment import run_real_pair

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--full-real',action='store_true',help='Use the full medium recording instead of 20 epochs')
args=parser.parse_args()
nb=nbformat.read('KFV-FGO-Colab.ipynb',as_version=4);nbformat.validate(nb)
assert len(nb.cells)==4
ns={}
for i,cell in enumerate(nb.cells[1:],1):
    with patch('IPython.display.display'):
        exec(compile(cell.source,f'colab_cell_{i}','exec'),ns)
    print('PASS widget module',i)
choices=[('EKF','riEKF'),('iEKF','Re-FGO'),('SW-FGO','SW-FGO')]
folder=Path('artifacts/colab-unified-smoke');folder.mkdir(parents=True,exist_ok=True)
for module,names in zip(MODULES,choices):
    for source in ['Simulation','Real GNSS/INS']:
        data,left,right=prepare_comparison(module,source,
            dict(method=names[0],iterations=2,kernel='huber'),
            dict(method=names[1],iterations=2,kernel='huber',window=3),
            simulation=dict(num_steps=30),match_filter_settings=module=='KFV vs FGO')
        if source=='Simulation': experiment=run_pair(data,left,right)
        else:
            if not args.full_real:
                data=copy(data);data.num_steps=20
                data.timestamps=data.timestamps[:20].copy()
            experiment=run_real_pair(data,left,right)
        preview=animate_pair(experiment).data
        assert '__PAYLOAD__' not in preview
        assert not experiment.get('groups')
        name=module.replace(' ','-').lower()+('-real' if source!='Simulation' else '-sim')
        (folder/(name+'.html')).write_text('<!doctype html><meta charset="utf-8">'+preview)
        print('PASS calculation',module,source,flush=True)
print('Six calculation paths and three notebook modules passed. Previews:',folder)
