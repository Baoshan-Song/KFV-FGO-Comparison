"""Reproducible independent KFV/Re-FGO equivalence checks, not accuracy tuning."""
import argparse
from datetime import datetime, timezone
import json
import platform
from pathlib import Path
import numpy as np
import scipy
from kfv_fgo import Config, FilterConfig, convert_kfv_to_fgo, generate_data, KfvEstimator, FgoEstimator


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path('artifacts/validation.json'))
    args=parser.parse_args()
    rows=[]
    for radius in (1000.,105.):
        for fraction in (0.,.2):
            data=generate_data(num_steps=100,emitter_radius=radius,gmm_weights=(1-fraction,fraction),seed=7)
            for mode in ('EKF','iEKF','rEKF','riEKF'):
                cfg=convert_kfv_to_fgo(Config(kfv=FilterConfig(mode=mode,max_iteration=10)))
                a,b=KfvEstimator(cfg,data).run(),FgoEstimator(cfg,data).run()
                state_error=float(np.max(np.abs(a['X_online']-b['X_online'])))
                covariance_error=float(np.max(np.abs(a['covariances']-b['covariances'])))
                np.testing.assert_allclose(a['X_online'],b['X_online'],atol=1e-8,rtol=0)
                np.testing.assert_allclose(a['covariances'],b['covariances'],atol=1e-10,rtol=1e-8)
                rows.append(dict(anchor_radius_m=radius,outlier_fraction=fraction,mode=mode,
                    max_iterations=cfg.fgo.max_iteration,range_sigma_m=.1,
                    max_state_difference=state_error,max_covariance_difference=covariance_error,
                    processed_epochs=99,step_tolerance_reached_epochs=sum(bool(d['converged']) for d in b['epoch_diagnostics'][1:])))
    report=dict(generated_utc=datetime.now(timezone.utc).isoformat(),python=platform.python_version(),
        numpy=np.__version__,scipy=scipy.__version__,platform=platform.platform(),seed=7,
        interpretation='Independent finite-iteration equivalence, not proof of convergence or reproduction of paper accuracy tables.',
        state_atol=1e-8,covariance_atol=1e-10,covariance_rtol=1e-8,cases=rows)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(f'{len(rows)} cases passed; max state difference {max(r["max_state_difference"] for r in rows):.3g}; max covariance difference {max(r["max_covariance_difference"] for r in rows):.3g}')
    print(args.output)

if __name__=='__main__':main()
