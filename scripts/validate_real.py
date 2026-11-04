"""Full-recording, independent real-data equivalence report."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import json
import numpy as np
from config.config import real_config,convert_kfv_to_fgo
from core.estimator import KfvEstimator,FgoEstimator
from kfv_fgo.data.gnss_ins import GnssImuDataset
cfg=real_config();data=GnssImuDataset(cfg)
report={'dataset':'urban_nav_deep','epochs':data.num_steps,
    'state_order':'ECEF position (m), ECEF velocity (m/s), body accel bias (m/s^2), receiver clock bias (m)',
    'solvers':'Independent gain update and graph QR; no filter dispatch from FGO',
    'modes':{}}
for mode in ('EKF','iEKF','rEKF','riEKF'):
    cfg.kfv.mode=mode;cfg=convert_kfv_to_fgo(cfg)
    a=KfvEstimator(cfg,data).run();b=FgoEstimator(cfg,data).run()
    diff=a['X_online']-b['X_online']
    report['modes'][mode]={
        'max_position_difference_m':float(np.max(np.linalg.norm(diff[:3],axis=0))),
        'max_abs_difference_per_state':np.max(abs(diff),axis=1).tolist(),
        'max_abs_covariance_difference':float(np.max(abs(a['covariances']-b['covariances']))),
        'state_check_atol':2e-6,'covariance_check_atol':2e-6}
    np.testing.assert_allclose(a['X_online'],b['X_online'],rtol=0,atol=2e-6)
    np.testing.assert_allclose(a['covariances'],b['covariances'],rtol=1e-7,atol=2e-6)
Path('docs/real_validation.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({m:r['max_position_difference_m'] for m,r in report['modes'].items()},indent=2))
