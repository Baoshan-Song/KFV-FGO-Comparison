"""Compute bundled Colab presets with the same independent production solvers."""
import argparse
import json
from pathlib import Path
import sys
from time import perf_counter

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from colab_ui import ComparisonUI
from colab_comparison import MODULES, prepare_comparison
from core.estimator import KfvEstimator, FgoEstimator
from kfv_fgo.data.gnss_ins import GnssImuDataset
from real_experiment import assemble_real_experiment
from run_cache import BUNDLED, method_key, save_result_cache, solver_fingerprint


def presets(dataset):
    unique={}
    for module in MODULES:
        ui=ComparisonUI(module)
        choices=([None,'w20'] if module=='FGO vs FGO' else
                 [None,'huber','cauchy','tukey'] if module=='KFV vs FGO' else [None])
        for kernel in choices:
            if kernel=='w20':ui.right.window.value=20
            elif kernel:
                ui.left.method.value='riEKF';ui.left.kernel.value=kernel
            _,left,right=prepare_comparison(module,'Real GNSS/INS',ui.left.values(),ui.right.values(),
                real=dict(dataset=dataset),match_filter_settings=ui.match_filter_settings,load_observations=False)
            for method in (left,right):unique[method_key(method)]=method
    return list(unique.values())


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--datasets',nargs='+',default=['medium','deep','harsh'])
    parser.add_argument('--output',type=Path,default=BUNDLED)
    parser.add_argument('--refresh',action='store_true',help='Recompute existing exact matches too')
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    manifest=args.output/'manifest.json'
    previous={r['file']:r for r in json.loads(manifest.read_text()).get('results',[])} if manifest.exists() else {}
    report=[]
    for dataset in args.datasets:
        methods=presets(dataset);data=GnssImuDataset.load(methods[0]['config'],progress=lambda x:print(dataset,x,flush=True))
        for method in methods:
            filename=method_key(method)+'.json.gz'
            if not args.refresh and filename in previous and (args.output/filename).is_file():
                report.append(previous[filename]);continue
            cfg=getattr(method['config'],method['kind']);started=perf_counter()
            result=(KfvEstimator if method['kind']=='kfv' else FgoEstimator)(method['config'],data).run()
            experiment=assemble_real_experiment(data,[method],[result])
            paths=save_result_cache(experiment,args.output)
            entry=dict(dataset=dataset,name=method['name'],kind=method['kind'],kernel=cfg.robust_kernel,
                iterations=cfg.max_iteration,window=cfg.window_size,delta=cfg.robust_delta,
                continuation=cfg.robust_continuation,epochs=data.num_steps-1,
                solver_seconds=result['runtime_seconds'],total_seconds=perf_counter()-started,
                bytes=paths[0].stat().st_size,file=paths[0].name)
            report.append(entry);print(json.dumps(entry),flush=True)
    (args.output/'manifest.json').write_text(json.dumps(dict(schema=1,solver=solver_fingerprint(),
        description='Independent full-recording runs; exact input/config/source cache matching.',results=report),indent=2)+'\n')


if __name__=='__main__':main()
