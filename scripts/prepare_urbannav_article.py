"""Prepare medium/harsh recordings from their public UrbanNav CSV/RINEX files.

The compact IMU archive preserves exactly the float64 columns read by ImuData.
No resampling, truth-derived attitude, time adjustment or coordinate fitting.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import numpy as np

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--navigation',type=Path,required=True)
    parser.add_argument('--scene',choices=['medium','harsh'],required=True)
    args=parser.parse_args()
    dest=Path(__file__).resolve().parents[1]/'kfv_fgo/data'/('urban_nav_'+args.scene)
    dest.mkdir(parents=True,exist_ok=True)
    files={'f9p_navi.obs':next(args.source.glob('*.ublox.f9p.obs')),
           'brdm.rnx':args.navigation,'gt.txt':args.source/'reference.txt'}
    for name,path in files.items():shutil.copyfile(path,dest/name)
    imu=args.source/'imu.csv'
    incomplete_tail=[]
    def complete_rows():
        with imu.open() as f:
            next(f)
            for number,line in enumerate(f,2):
                if len(line.split(',')) < 32:
                    if f.read():raise ValueError(f'Incomplete IMU row {number} is not the final row')
                    incomplete_tail.append(number)
                    break
                yield line
    v=np.loadtxt(complete_rows(),delimiter=',',usecols=(0,4,5,6,7,29,30,31))
    np.savez_compressed(dest/'xsens_imu.npz',time_s=v[:,0]*1e-9,quat_xyzw=v[:,1:5],acc_mps2=v[:,5:8])
    files['xsens_imu.npz']=imu
    metadata=dict(sequence=args.source.name,source='https://github.com/IPNL-POLYU/UrbanNavDataset',
        conversion='IMU: retain all rows and the eight float64 columns consumed by the original CSV loader, converting nanoseconds to seconds. All other files are byte-identical copies.',
        source_sha256={name:sha(path) for name,path in files.items()},
        packaged_sha256={name:sha(dest/name) for name in files},imu_rows=len(v),
        incomplete_final_rows_omitted=incomplete_tail)
    (dest/'provenance.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print(dest, len(v),'IMU rows',flush=True)
