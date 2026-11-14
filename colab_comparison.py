"""Shared configuration and execution for the three Colab comparison modules."""
from copy import deepcopy
from pathlib import Path

import numpy as np

from config.config import convert_kfv_to_fgo, real_config
from experiment_common import make_config, make_data, print_statistics, run_pair


MODULES = {'KFV vs KFV': ('kfv', 'kfv'),
           'KFV vs FGO': ('kfv', 'fgo'),
           'FGO vs FGO': ('fgo', 'fgo')}
FILTERS = ('EKF', 'iEKF', 'rEKF', 'riEKF')
GRAPHS = ('Re-FGO', 'SW-FGO')


def _configured_method(kind, settings, base):
    cfg = deepcopy(base)
    name = settings['method']
    allowed = FILTERS if kind == 'kfv' else GRAPHS
    if name not in allowed:
        raise ValueError(f'{name} is not a {kind.upper()} method')
    iterations = int(settings.get('iterations', 5))
    window = int(settings.get('window', 5))
    delta = float(settings.get('delta', 1.))
    kernel = settings.get('kernel', 'none')
    if iterations < 1 or window < 1:
        raise ValueError('Iterations and window size must be positive')
    if not np.isfinite(delta) or delta <= 0:
        raise ValueError('Robust scale must be finite and positive')
    if kernel not in ('none', 'huber', 'cauchy', 'tukey'):
        raise ValueError(f'Unsupported robust kernel: {kernel}')
    getattr(cfg, kind).robust_continuation = bool(settings.get('continuation', False))
    if kind == 'kfv':
        cfg.kfv.mode = name
        cfg.kfv.max_iteration = iterations
        cfg.kfv.robust_kernel = kernel
        cfg.kfv.robust_delta = delta
    else:
        cfg.fgo.max_iteration = iterations
        cfg.fgo.robust_kernel = kernel
        cfg.fgo.robust_delta = delta
        cfg.fgo.imitate_kfv = name == 'Re-FGO'
        cfg.fgo.window_size = 1 if cfg.fgo.imitate_kfv else window
        cfg.fgo.autodiff = bool(settings.get('numerical_jacobian', False)) if cfg.data_mode == 'sim' else False
        if name == 'SW-FGO':
            name = f'SW-FGO (W={window})'
    return dict(name=name, kind=kind, config=cfg)


def prepare_comparison(module, data_source, left, right, *, simulation=None,
                       real=None, match_filter_settings=False, load_observations=True,
                       progress=None, stop=None):
    """Build two independent configurations using the same observations and prior."""
    if module not in MODULES:
        raise ValueError(f'Unknown comparison module: {module}')
    if data_source not in ('Simulation', 'Real GNSS/INS'):
        raise ValueError(f'Unknown data source: {data_source}')
    if match_filter_settings and (module != 'KFV vs FGO' or right['method'] != 'Re-FGO'):
        raise ValueError('Match filter settings requires KFV vs FGO with Re-FGO on the right')
    if data_source == 'Simulation':
        options = dict(simulation or {})
        err_x = options.pop('initial_error_x', 100.)
        err_y = options.pop('initial_error_y', -100.)
        sigma = float(options.get('white_sigma', .1))
        if not np.isfinite(sigma) or sigma <= 0:
            raise ValueError('Nominal range noise must be finite and positive')
        base = make_config(err_x=err_x, err_y=err_y)
        base.kfv.r = base.fgo.r = sigma ** 2
        if progress: progress('Generating simulated observations')
        data = make_data(**options) if load_observations else None
    else:
        from kfv_fgo.data.gnss_ins import GnssImuDataset
        options = dict(real or {})
        base = real_config()
        directory = options.get('data_directory', '').strip()
        base.data_path = str(Path(directory).expanduser().resolve()) if directory else 'urban_nav_' + options.get('dataset', 'medium')
        base.leap_seconds = int(options.get('gps_minus_utc_s', 18))
        elevation = float(options.get('minimum_elevation_deg', -90.))
        snr = float(options.get('minimum_snr_dbhz', 0.))
        base.gnss['minimum_elevation_deg'] = None if elevation <= -90 else elevation
        base.gnss['minimum_snr'] = None if snr <= 0 else snr
        std = np.array([options.get('position_process_std_m', .3)] * 3 +
                       [options.get('velocity_process_std_mps', .15)] * 3 +
                       [options.get('bias_process_std_mps2', .01)] * 3 +
                       [options.get('clock_process_std_m', 100.)], dtype=float)
        if not np.all(np.isfinite(std)) or np.any(std <= 0):
            raise ValueError('Process standard deviations must be finite and positive')
        base.kfv.q = np.diag(std ** 2)
        base = convert_kfv_to_fgo(base)
        data = GnssImuDataset.load(base, progress=progress, stop=stop) if load_observations else None
    first, second = [_configured_method(kind, settings, base)
                     for kind, settings in zip(MODULES[module], (left, right))]
    if match_filter_settings:
        second['config'] = convert_kfv_to_fgo(first['config'])
        second['name'] = f"Re-FGO (matched {first['name']})"
    return data, first, second


def run_comparison(module, data_source, left, right, *, output_directory='artifacts/colab', **options):
    data, first, second = prepare_comparison(module, data_source, left, right, **options)
    if data_source == 'Simulation':
        experiment = run_pair(data, first, second)
        print_statistics({f'{i + 1}. {item["name"]}': item['metrics']
                          for i, item in enumerate(experiment['estimators'])})
    else:
        from real_experiment import run_real_pair, save_real_experiment
        print(f'{data.num_steps} GNSS epochs; UTC duration {np.ptp(data.timestamps):.1f} s')
        experiment = run_real_pair(data, first, second)
        for item in experiment['estimators']:
            report = item['accuracy']
            if report:
                print(f"{item['name']}: H {report['horizontal']['rmse']:.3f} m; "
                      f"Up {report['ENU']['up']['rmse']:.3f} m; "
                      f"3D {report['three_dimensional']['rmse']:.3f} m")
        # Preserve each run rather than silently overwriting another comparison.
        from datetime import datetime
        folder = Path(output_directory) / module.lower().replace(' ', '-') / datetime.now().strftime('%Y%m%d-%H%M%S-%f')
        save_real_experiment(experiment, folder)
        print('Saved viewer, configurations, states and diagnostics to', folder.resolve())
    return experiment
