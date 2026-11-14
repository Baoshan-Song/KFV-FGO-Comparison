"""Responsive partial-result viewer and epoch-interleaved comparison runner."""
from copy import deepcopy
from html import escape
from threading import RLock
from time import monotonic
from urllib.parse import quote

import numpy as np
import ipywidgets as W

from core.estimator import KfvEstimator, FgoEstimator
from experiment_common import metrics, measurement_residuals


def run_interleaved(data, methods, receive, stop):
    """Yield real solver updates, alternating methods at each observation epoch.

    UI work runs outside measured epoch updates. The numerical paths are the
    same generators consumed by ordinary Estimator.run().
    """
    streams=[]
    for item in methods:
        cfg=deepcopy(item['config'])
        if isinstance(data,dict):
            for c in (cfg.kfv,cfg.fgo):
                c.dt=float(np.asarray(data['dt']).squeeze())
                c.omega=float(np.asarray(data['omega']).squeeze())
        streams.append((KfvEstimator if item['kind']=='kfv' else FgoEstimator)(cfg,data).iter_epochs())
    results=[None]*len(streams)
    try:
        while any(result is None for result in results):
            for i,stream in enumerate(streams):
                if stop.is_set():return None
                if results[i] is not None:continue
                try: receive(i,next(stream))
                except StopIteration as completed:results[i]=completed.value
    finally:
        for stream in streams:stream.close()
    if not isinstance(data,dict):
        from real_experiment import assemble_real_experiment
        return assemble_real_experiment(data,methods,results)
    return dict(data=data,estimators=[{**item,'result':result,
        'residuals':measurement_residuals(result,data),
        'metrics':metrics(result,data['true_positions'],result['runtime_seconds']*1000)}
        for item,result in zip(methods,results)])


class LiveComparison:
    def __init__(self,data,methods):
        self.data,self.names=data,[m['name'] for m in methods]
        self.real=not isinstance(data,dict)
        self.count=data.num_steps if self.real else data['num_steps']
        self.time=(data.timestamps-data.timestamps[0] if self.real else
                   np.arange(self.count)*float(np.asarray(data['dt']).squeeze()))
        self.gap_limit=2*float(np.median(np.diff(self.time)))
        self.dim=3 if self.real else 2
        self.position=np.full((2,self.count,self.dim),np.nan)
        self.error=np.full_like(self.position,np.nan)
        self.norms=np.full((2,2,self.count),np.nan)
        self.ready=[0,0];self.lock=RLock();self._internal=False;self._last_render=0.
        self.reference=np.full((self.count,self.dim),np.nan)
        if self.real:
            from kfv_fgo.model.gnss_ins import local_to_ecef
            from kfv_fgo.core.evaluation import load_ground_truth
            self.origin=data.initial_state[:3].copy();self.rotation=local_to_ecef(self.origin)
            self.truth_ecef=np.full((self.count,3),np.nan)
            path=data.path/'gt.txt'
            if path.exists():
                t,xyz=load_ground_truth(path)
                valid=(data.timestamps>=t[0])&(data.timestamps<=t[-1]);valid[0]=False
                if valid.any():
                    self.truth_ecef[valid]=np.column_stack([np.interp(data.timestamps[valid],t,xyz[:,j]) for j in range(3)])
                    self.error_rotation=local_to_ecef(self.truth_ecef[valid][0])
                    self.reference=(self.truth_ecef-self.origin)@self.rotation
        else:self.reference=data['true_positions'].T.copy()
        self.has_truth=np.isfinite(self.reference[1:]).any()
        self.view=W.ToggleButtons(options=['Trajectory','Position error','Residual norm'])
        self.component=W.Dropdown(options=['Horizontal','Up','3D'] if self.real else ['Horizontal'],layout=W.Layout(width='150px'))
        self.residual=W.Dropdown(options=['Post-fit','Innovation'],layout=W.Layout(width='150px'))
        self.follow=W.Checkbox(value=True,description='Follow latest',indent=False)
        self.timeline=W.IntSlider(value=0,min=0,max=1,description='Time',continuous_update=True,
                                   readout=False,layout=W.Layout(width='70%'))
        self.detail=W.IntSlider(value=30,min=5,max=120,description='Detail (s)',continuous_update=False,
                                layout=W.Layout(width='280px'))
        self.clock=W.HTML();self.figure=W.HTML();self.download=W.HTML()
        self.progress=W.IntProgress(min=0,max=2*(self.count-1),description='Solved',layout=W.Layout(width='95%'))
        self.note=W.HTML()
        self.widget=W.VBox([self.progress,self.view,W.HBox([self.component,self.residual,self.follow]),
            W.HBox([self.timeline,self.clock]),self.detail,self.note,self.figure,self.download])
        self.timeline.observe(self._seek,names='value')
        for widget in (self.view,self.component,self.residual,self.detail,self.follow):widget.observe(self._control,names='value')
        self.render()

    def _seek(self,_):
        if self._internal:return
        self.follow.value=False
        self.render()

    def _control(self,_):
        if self._internal:return
        self.render()

    def receive(self,method,event):
        with self.lock:
            k=event['epoch'];x=event['state'];d=event['diagnostic']
            if self.real:
                self.position[method,k]=(x[:3]-self.origin)@self.rotation
                if np.isfinite(self.truth_ecef[k]).all():
                    self.error[method,k]=(x[:3]-self.truth_ecef[k])@self.error_rotation
                self.norms[method,0,k]=np.linalg.norm(d['postfit_residual'])
                self.norms[method,1,k]=np.linalg.norm(d['prefit_residual'])
            else:
                self.position[method,k]=x[:2]
                self.error[method,k]=x[:2]-self.reference[k]
                z=self.data['toa_measurements'][:,k];emitters=self.data['emitter_positions']
                for j,state in enumerate((x,event['predicted'])):
                    self.norms[method,j,k]=np.linalg.norm(z-np.linalg.norm(state[:2,None]-emitters,axis=0))
            self.ready[method]=k
        if monotonic()-self._last_render>.3:self.render()

    def finish(self,stopped=False):
        self.progress.bar_style='warning' if stopped else 'success'
        self.render()

    def render(self):
        with self.lock:
            self._last_render=monotonic();available=min(self.ready)
            self._internal=True
            try:
                self.progress.value=sum(self.ready)
                self.timeline.max=max(1,available)
                if self.follow.value:self.timeline.value=available
                k=min(self.timeline.value,available)
                self.clock.value=f'<b>{self.time[k]:.1f} s</b> · epoch {k}'
                self.component.layout.display='' if self.view.value=='Position error' else 'none'
                self.residual.layout.display='' if self.view.value=='Residual norm' else 'none'
            finally:self._internal=False
            self.note.value=(f'<small>Computed through {self.time[available]:.1f} s. Drag the time slider to inspect completed epochs; '
                             'enable Follow latest to resume tracking. Grey bands mark recorded GNSS gaps.</small>')
            if available<1:
                self.figure.value='<p>Waiting for the first pair of updates…</p>';return
            all_indices=np.arange(1,available+1)
            detail=all_indices[(self.time[all_indices]>=self.time[k]-self.detail.value)&(all_indices<=k)]
            if not len(detail):detail=np.array([1])
            kind=self.view.value
            if kind=='Trajectory':
                values=self.position[:,:,:2];reference=self.reference[:,:2]
                xlabel,ylabel=('East (m)','North (m)') if self.real else ('X (m)','Y (m)')
            else:
                if kind=='Residual norm':
                    scalar=self.norms[:,0 if self.residual.value=='Post-fit' else 1]
                    ylabel=f'{self.residual.value} residual norm (m)'
                elif self.component.value=='Up':scalar=self.error[:,:,2];ylabel='Up error (m)'
                elif self.component.value=='3D':scalar=np.linalg.norm(self.error,axis=2);ylabel='3D error (m)'
                else:scalar=np.linalg.norm(self.error[:,:,:2],axis=2);ylabel='Horizontal error (m)'
                values=np.stack([np.broadcast_to(self.time,scalar.shape),scalar],axis=-1)
                reference=None;xlabel='Elapsed time (s)'
            svg=['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1080 410" style="width:100%;background:white;font-family:Arial;font-size:13px">']
            colors=['#1857d3','#e3272b']
            for i,name in enumerate(self.names):svg.append(f'<text x="{90+i*380}" y="23" fill="{colors[i]}">{escape(name)}</text>')
            if reference is not None:svg.append('<text x="850" y="23" fill="#555">Reference</text>')
            for panel,indices in enumerate((all_indices,detail)):
                left=75+panel*535;top=65;width=440;height=275
                sets=[v[indices] for v in values]
                if reference is not None:sets.insert(0,reference[indices])
                points=np.concatenate(sets);valid=points[np.isfinite(points).all(axis=1)]
                if not len(valid):
                    svg.append(f'<text x="{left}" y="160">No overlapping reference data yet</text>');continue
                lo=valid.min(axis=0);hi=valid.max(axis=0);span=np.maximum(hi-lo,.1)
                if kind=='Trajectory':span=np.array([width,height])*max(span/[width,height])
                mid=(lo+hi)/2;lo=mid-span*.55;hi=mid+span*.55
                def xy(p):return (left+(p[0]-lo[0])/(hi[0]-lo[0])*width,top+height-(p[1]-lo[1])/(hi[1]-lo[1])*height)
                svg.append(f'<text x="{left}" y="49">{ "Computed interval" if panel==0 else "Selected time window"}</text>')
                for frac in np.linspace(0,1,5):
                    X=left+frac*width;Y=top+height-frac*height
                    svg.append(f'<path d="M{X},{top}v{height} M{left},{Y}h{width}" stroke="#e4e9ef" fill="none"/>')
                    svg.append(f'<text x="{X}" y="{top+height+21}" text-anchor="middle">{lo[0]+frac*(hi[0]-lo[0]):.1f}</text>')
                    svg.append(f'<text x="{left-8}" y="{Y+4}" text-anchor="end">{lo[1]+frac*(hi[1]-lo[1]):.1f}</text>')
                for i,curve in enumerate(sets):
                    color='#555' if reference is not None and i==0 else colors[i-(reference is not None)]
                    segments=[];active=False;previous=None
                    for epoch,p in zip(indices,curve):
                        if not np.isfinite(p).all():active=False;continue
                        if previous is not None and self.time[epoch]-self.time[previous]>self.gap_limit:
                            active=False
                            if kind!='Trajectory':
                                x1,_=xy([self.time[previous],0]);x2,_=xy([self.time[epoch],0])
                                svg.append(f'<rect x="{x1}" y="{top}" width="{x2-x1}" height="{height}" fill="#aaa" opacity=".1"/>')
                        X,Y=xy(p);segments.append(f'{"L" if active else "M"}{X:.2f},{Y:.2f}');active=True;previous=epoch
                    svg.append(f'<path d="{" ".join(segments)}" stroke="{color}" stroke-width="2" fill="none"/>')
                if kind!='Trajectory':
                    X,_=xy([self.time[k],0]);svg.append(f'<path d="M{X},{top}v{height}" stroke="#555" stroke-dasharray="4 4"/>')
                else:
                    for i in range(2):
                        if np.isfinite(values[i,k]).all():
                            X,Y=xy(values[i,k]);svg.append(f'<circle cx="{X}" cy="{Y}" r="5" fill="{colors[i]}"/>')
                svg.append(f'<text x="{left+width/2}" y="392" text-anchor="middle">{xlabel}</text>')
                svg.append(f'<text transform="translate({left-55},{top+height/2}) rotate(-90)" text-anchor="middle">{ylabel}</text>')
            svg.append('</svg>');self.svg=''.join(svg);self.figure.value=self.svg
            self.download.value=f'<a download="comparison-live.svg" href="data:image/svg+xml;charset=utf-8,{quote(self.svg)}">Save visible SVG</a>'
