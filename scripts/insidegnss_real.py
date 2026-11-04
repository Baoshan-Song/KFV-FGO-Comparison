"""Produce the article's real-data runs with the SAME Colab presets and viewer.

Run python scripts/insidegnss_real.py. Open each viewer.html and use Save SVG.
Figures therefore share exactly the browser renderer used in Colab.
"""
import argparse
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import json
import platform
import numpy as np
from real_experiment import preset_pair,run_real_pair,save_real_experiment,clean
from config.config import real_config
from kfv_fgo.data.gnss_ins import GnssImuDataset
from kfv_fgo.core.evaluation import compare_results


def run(output):
    output=Path(output)
    data=GnssImuDataset(real_config())
    summary={'epochs':data.num_steps,'platform':platform.platform(),
             'python':platform.python_version(),'numpy':np.__version__,
             'evaluated_epochs':None,'presets':{}}
    for preset in ('Matched EKF','Matched riEKF','EKF vs riEKF','Re-FGO vs SW-FGO'):
        first,second=preset_pair(preset)
        experiment=run_real_pair(data,first,second)
        directory=output/preset.lower().replace(' ','-')
        save_real_experiment(experiment,directory)
        rows={}
        for item in experiment['estimators']:
            r=item['result'];a=item['accuracy']
            rows[item['name']]={'horizontal_rmse_m':a['horizontal']['rmse'],
                'up_rmse_m':a['ENU']['up']['rmse'],
                'spatial_rmse_m':a['three_dimensional']['rmse'],
                'ms_per_update':r['runtime_seconds']*1000/(data.num_steps-1)}
            summary['evaluated_epochs']=a['evaluated_epochs']
        a,b=[i['result'] for i in experiment['estimators']]
        summary['presets'][preset]={'methods':rows,'difference':compare_results(
            {**a,'X':a['X_online']},{**b,'X':b['X_online']})}
        print(preset,rows,flush=True)
    output.mkdir(exist_ok=True,parents=True)
    (output/'article-summary.json').write_text(json.dumps(clean(summary),indent=2,allow_nan=False)+'\n')

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path('artifacts/insidegnss-real'))
    run(parser.parse_args().output)
