"""Six-method Inside GNSS comparison, identical to the Colab article cell."""
import argparse
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from config.config import real_config
from kfv_fgo.data.gnss_ins import GnssImuDataset
from real_experiment import run_article_comparison, save_real_experiment, clean

def cpu_name():
    if sys.platform=='darwin':
        try:
            return subprocess.check_output(['sysctl','-n','machdep.cpu.brand_string'],stderr=subprocess.DEVNULL,text=True).strip()
        except (OSError,subprocess.CalledProcessError):
            pass
    return platform.processor() or platform.machine()

if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scene',choices=['medium','deep','harsh'],default='medium')
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    if args.output is None:args.output=Path('artifacts/insidegnss-real')/('article-'+args.scene+'-w20')
    cfg=real_config();cfg.data_path='urban_nav_'+args.scene
    data=GnssImuDataset(cfg)
    experiment=run_article_comparison(data,cfg)
    save_real_experiment(experiment,args.output)
    rows=[]
    for item in experiment['estimators']:
        a,r=item['accuracy'],item['result'];t=r['epoch_runtime_ms'][1:]
        rows.append(dict(name=item['name'],H=a['horizontal']['rmse'],Up=a['ENU']['up']['rmse'],
            D3=a['three_dimensional']['rmse'],epochs=a['evaluated_epochs'],
            median_ms=float(np.median(t)),q1_ms=float(np.quantile(t,.25)),q3_ms=float(np.quantile(t,.75)),
            mean_ms=float(np.mean(t)),min_ms=float(t.min()),max_ms=float(t.max())))
    diffs={}
    for label,offset in [('nonrobust',0),('robust',3)]:
        a,b=[experiment['estimators'][offset+j]['result']['X_online'][:3] for j in (0,1)]
        diffs[label]=float(np.max(np.linalg.norm(a-b,axis=0)))
    provenance=dict(scene=args.scene,rows=rows,matched_max_3d_m=diffs,platform=platform.platform(),
        python=platform.python_version(),numpy=np.__version__,
        cpu=cpu_name(),
        threads={k:os.getenv(k) for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS')},
        timing=f'Serial full recordings, 25-epoch warmup discarded; {data.num_steps-1} full-update timings per method, including propagation, diagnostics and marginalization; excluding load, warmup, evaluation and render.',
        source_revision=subprocess.getoutput('git rev-parse HEAD'),
        source_diff=subprocess.getoutput('git diff -- core/estimator/real_data.py real_experiment.py'))
    (args.output/'article-summary.json').write_text(json.dumps(clean(provenance),indent=2))
    print(json.dumps(clean(provenance),indent=2))
