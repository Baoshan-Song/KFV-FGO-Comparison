"""Content-addressed results for the Colab viewer. No pickle or estimator shortcut.

Only an exact match of inputs, numerical configuration and solver source can
reuse a result. A changed parameter or input starts an independent fresh solve.
"""
from copy import deepcopy
from functools import lru_cache
import gzip
import hashlib
import json
from pathlib import Path
import platform
from datetime import datetime, timezone

from config.config import _config_to_dict, resolve_data_path

ROOT=Path(__file__).resolve().parent
BUNDLED=ROOT/'kfv_fgo'/'precomputed'
LOCAL=ROOT/'artifacts'/'comparison-cache'
SCHEMA=1


@lru_cache(maxsize=128)
def _file_hash(path,size,mtime_ns):
    h=hashlib.sha256()
    with open(path,'rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def _hash_file(path):
    stat=path.stat()
    return _file_hash(str(path),stat.st_size,stat.st_mtime_ns)


def _digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def input_fingerprint(config):
    path=resolve_data_path(config)
    names=['f9p_navi.obs','brdm.rnx','xsens_imu.csv' if (path/'xsens_imu.csv').is_file() else 'xsens_imu.npz']
    if (path/'gt.txt').is_file():names.append('gt.txt')
    return _digest({name:_hash_file(path/name) for name in names})


def solver_fingerprint():
    paths=list((ROOT/'core').rglob('*.py'))+list((ROOT/'config').glob('*.py'))
    paths += [ROOT/p for p in ('kfv_fgo/model/gnss_ins.py','kfv_fgo/data/gnss_ins.py',
        'kfv_fgo/data/rtklib.py','kfv_fgo/core/evaluation.py','data/circle_eval.py','real_experiment.py')]
    return _digest({str(p.relative_to(ROOT)):_hash_file(p) for p in sorted(paths)})


def method_key(method,inputs=None):
    cfg=method['config'];record=_config_to_dict(cfg)
    record.pop('FGO' if method['kind']=='kfv' else 'KFV')
    # The content hash defines the recording, independent of its checkout path.
    record['data']['path']='content-addressed'
    return _digest(dict(schema=SCHEMA,solver=solver_fingerprint(),kind=method['kind'],config=record,
                        inputs=inputs if inputs is not None else input_fingerprint(cfg)))


PER_METHOD=('names','kinds','trajectories','diagnostics','metrics','errors','satellite_counts','bias','epoch_runtime_ms')


def save_result_cache(experiment,directory=LOCAL):
    """Store each independently computed method for reuse across UI modules."""
    from real_experiment import real_payload
    if isinstance(experiment['data'],dict):return
    payload=real_payload(experiment);directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    inputs=input_fingerprint(experiment['estimators'][0]['config'])
    files=[]
    for i,method in enumerate(experiment['estimators']):
        key=method_key(method,inputs)
        single={k:v for k,v in payload.items() if k not in PER_METHOD and k not in ('residuals','metadata','groups')}
        single.update({k:[payload[k][i]] for k in PER_METHOD})
        single['residuals']={k:[v[i]] for k,v in payload['residuals'].items()}
        single['metadata']={**payload['metadata'],'configs':[payload['metadata']['configs'][i]]}
        single['groups']=[]
        record=dict(schema=SCHEMA,key=key,payload=single,computed_at=datetime.now(timezone.utc).isoformat(),
                    machine=f'{platform.system()} {platform.machine()}',solver=solver_fingerprint())
        target=directory/(key+'.json.gz');temp=target.with_suffix('.tmp')
        temp.write_bytes(gzip.compress(json.dumps(record,separators=(',',':'),allow_nan=False).encode(),mtime=0))
        temp.replace(target);files.append(target)
    return files


def load_result_cache(methods,directories=None):
    directories=directories or (LOCAL,BUNDLED)
    inputs=input_fingerprint(methods[0]['config']);records=[]
    for method in methods:
        key=method_key(method,inputs);record=None
        for directory in directories:
            path=Path(directory)/(key+'.json.gz')
            if not path.is_file():continue
            try:
                candidate=json.loads(gzip.decompress(path.read_bytes()))
                if candidate.get('schema')==SCHEMA and candidate.get('key')==key:
                    record=candidate;break
            except (OSError,ValueError,EOFError):continue
        if record is None:return None
        records.append(record)
    first=records[0]['payload'];payload=deepcopy(first)
    # Never splice trajectories sampled on different time grids.
    if any(r['payload']['time']!=first['time'] for r in records):return None
    for field in PER_METHOD:payload[field]=[r['payload'][field][0] for r in records]
    payload['names']=[m['name'] for m in methods]
    payload['residuals']={key:[r['payload']['residuals'][key][0] for r in records] for key in first['residuals']}
    payload['metadata']['configs']=[r['payload']['metadata']['configs'][0] for r in records]
    payload['metadata']['cache']=dict(reused=True,computed_at=[r['computed_at'] for r in records],
        machines=[r['machine'] for r in records],timing_note='Stored solve timings, not this cache load.')
    return payload
