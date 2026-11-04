"""Real GNSS/INS experiments and the payload used by the Colab SVG viewer."""
from copy import deepcopy
from pathlib import Path
import json
import numpy as np
from config.config import real_config, convert_kfv_to_fgo
from core.estimator import KfvEstimator, FgoEstimator
from kfv_fgo.data.gnss_ins import GnssImuDataset
from kfv_fgo.core.evaluation import load_ground_truth, evaluate_ground_truth
from kfv_fgo.model.gnss_ins import local_to_ecef

PRESETS = ('Matched EKF', 'Matched riEKF', 'EKF vs riEKF', 'Re-FGO vs SW-FGO', 'Custom')


def method(name, base, iterations=5, kernel='huber', delta=1., window=5, recursive_mode='riEKF'):
    cfg = deepcopy(base)
    mode = recursive_mode if name in ('Re-FGO', 'SW-FGO') else name
    if mode not in ('EKF','iEKF','rEKF','riEKF'):
        raise ValueError(f'Unsupported method {name}')
    cfg.kfv.mode, cfg.kfv.max_iteration = mode, int(iterations)
    cfg.kfv.robust_kernel, cfg.kfv.robust_delta = kernel, float(delta)
    cfg = convert_kfv_to_fgo(cfg)
    graph = name in ('Re-FGO', 'SW-FGO')
    if name == 'SW-FGO':
        cfg.fgo.imitate_kfv = False
        cfg.fgo.window_size = int(window)
    return dict(name=(f'{name} ({mode})' if name=='Re-FGO' else
                      f'SW-FGO (W={window})' if name=='SW-FGO' else name),
                kind='fgo' if graph else 'kfv', config=cfg)


def preset_pair(preset, base=None, **options):
    base = real_config() if base is None else base
    choices = {'Matched EKF':('EKF','Re-FGO','EKF'),
               'Matched riEKF':('riEKF','Re-FGO','riEKF'),
               'EKF vs riEKF':('EKF','riEKF','riEKF'),
               'Re-FGO vs SW-FGO':('Re-FGO','SW-FGO','riEKF')}
    a,b,mode=choices[preset]
    options['recursive_mode']=mode
    return method(a,base,**options), method(b,base,**options)


def clean(value):
    if isinstance(value,np.ndarray): return clean(value.tolist())
    if isinstance(value,np.generic): return clean(value.item())
    if isinstance(value,dict): return {str(k):clean(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)): return [clean(v) for v in value]
    if isinstance(value,float): return value if np.isfinite(value) else None
    return value


def run_real_pair(data, first, second):
    truth_path=data.path/'gt.txt'
    truth=load_ground_truth(truth_path) if truth_path.exists() else None
    estimators=[]
    for item in (first,second):
        print('Running',item['name'],flush=True)
        result=(KfvEstimator if item['kind']=='kfv' else FgoEstimator)(item['config'],data).run()
        # All viewer metrics are causal, including SW-FGO.
        causal={**result,'X':result['X_online']}
        report, errors = evaluate_ground_truth(causal,truth) if truth else (None,None)
        estimators.append({**item,'result':result,'accuracy':report,'errors':errors})
    return {'data':data,'estimators':estimators}


def real_payload(experiment):
    data,items=experiment['data'],experiment['estimators']
    n=data.num_steps
    origin=data.initial_state[:3]
    rotation=local_to_ecef(origin)
    truth=np.full((3,n),np.nan)
    if items[0]['errors'] is not None:
        e=items[0]['errors']
        truth[:,e['valid']]=((e['truth_ecef']-origin)@rotation).T
    sat_ids=sorted({int(s) for ep in data.epochs for s in ep.sat_ids})
    rows={s:i for i,s in enumerate(sat_ids)}
    trajectories,diagnostics,metrics,errors,bias,sat_counts=[],[],[],[],[],[]
    residuals={key:[] for key in ('postfit','prefit')}
    for item in items:
        result=item['result']; x=result['X_online']
        trajectories.append(((x[:3].T-origin)@rotation).T)
        traces=deepcopy(result['epoch_diagnostics'])
        counts=np.full(n,np.nan)
        for k,record in enumerate(traces):
            if record is None:continue
            counts[k]=len(record['postfit_ids'])
            for step in record['iteration_trace']:
                step['state'][:3]=(step['state'][:3]-origin)@rotation
        diagnostics.append(traces)
        sat_counts.append(counts)
        bias.append(x[6:9])
        for mode in residuals:
            values=np.full((len(sat_ids),n),np.nan)
            for k,record in enumerate(result['epoch_diagnostics']):
                if record is None:continue
                for s,r in zip(record[mode+'_ids'],record[mode+'_residual']):
                    values[rows[int(s)],k]=r
            residuals[mode].append(values)
        error={key:np.full(n,np.nan) for key in ('horizontal','up','3d')}
        if item['accuracy']:
            e,report=item['errors'],item['accuracy']
            error['horizontal'][e['valid']]=e['horizontal']
            error['up'][e['valid']]=e['error_enu'][:,2]
            error['3d'][e['valid']]=e['three_dimensional']
            h=report['horizontal']
            metric=dict(rmse=h['rmse'],up_rmse=report['ENU']['up']['rmse'],
                spatial_rmse=report['three_dimensional']['rmse'],
                mae=h['mae'],cp95=h['absolute_error_95'],max_error=h['max_error'])
        else:
            metric=dict(rmse=None,up_rmse=None,spatial_rmse=None)
        metric.update(elapsed_ms=result['runtime_seconds']*1000,
            per_step_ms=result['runtime_seconds']*1000/(n-1),processed_epochs=n-1)
        metrics.append(metric);errors.append(error)
    return clean(dict(real=True,names=[i['name'] for i in items],
        kinds=[i['kind'] for i in items],frames=n,dt=float(np.median(np.diff(data.timestamps))),
        time=data.timestamps-data.timestamps[0],truth=truth,trajectories=trajectories,
        residuals=residuals,diagnostics=diagnostics,metrics=metrics,errors=errors,
        satellite_counts=sat_counts,satellite_ids=sat_ids,bias=bias,
        has_truth=bool(items[0]['accuracy']),
        metadata=dict(model='10-state pseudorange/IMU with supplied AHRS attitude',
            state_order='ECEF position, ECEF velocity, body accelerometer bias, clock bias',
            time_scale='UTC',dataset=data.path.name,load_seconds=data.load_seconds,
            configs=[i['result']['config'] for i in items])))


def save_real_experiment(experiment, output):
    from experiment_common import animate_pair
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    payload=real_payload(experiment)
    (output/'viewer.html').write_text('<!doctype html><meta charset="utf-8">'+animate_pair(experiment).data)
    (output/'viewer-data.json').write_text(json.dumps(payload,allow_nan=False))
    summary={}
    for index,item in enumerate(experiment['estimators']):
        result=item['result'];label=f'{index+1}-{item["name"]}'.replace(' ','_')
        np.savez_compressed(output/f'{label}.npz',**{k:result[k] for k in (
            'X','X_online','X_pred','covariances','timestamps')})
        np.savetxt(output/f'{label}.csv',np.column_stack((result['timestamps'],result['X_online'].T)),
            delimiter=',',header='utc_s,x_m,y_m,z_m,vx_mps,vy_mps,vz_mps,bax_mps2,bay_mps2,baz_mps2,clock_m',comments='')
        summary[label]=dict(accuracy=item['accuracy'],runtime_seconds=result['runtime_seconds'],
            config=result['config'],solver=result['solver'])
    (output/'summary.json').write_text(json.dumps(clean(summary),indent=2,allow_nan=False))
    return output
