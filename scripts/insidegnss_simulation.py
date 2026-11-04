"""Export reproducible simulation figures in the Inside GNSS draft's style."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, ConnectionPatch
from matplotlib.ticker import MaxNLocator
from config.config import Config, FilterConfig, FgoConfig, convert_kfv_to_fgo
from core.estimator import KfvEstimator, FgoEstimator
from data.circle_eval import generate_data
from experiment_common import measurement_residuals

BLUE,RED='#1857d3','#e3272b'
STYLE={
 'EKF':('#1857d3','x'), 'rEKF':('#ce00ce','x'),
 'iEKF':('#b47700','x'), 'riEKF':('#16803a','x'),
 'Re-FGO':(RED,'o'), 'SW-FGO, W=5':('#00aabd','o')}


def appearance():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':13,
        'axes.labelsize':14,'axes.titlesize':14,'legend.fontsize':11.5,
        'lines.linewidth':1.9,'lines.markersize':5,'pdf.fonttype':42,
        'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,
        'savefig.facecolor':'white','figure.facecolor':'white'})


def axes(ax,xlabel,ylabel):
    ax.set(xlabel=xlabel,ylabel=ylabel)
    ax.grid(True,color='#dedede',linewidth=.6)
    ax.xaxis.set_major_locator(MaxNLocator(4)); ax.yaxis.set_major_locator(MaxNLocator(5))


def draw(ax,x,y,name,every=1):
    color,mark=STYLE[name]
    return ax.plot(x,y,color=color,marker=mark,markevery=every,
        markerfacecolor='none' if mark=='o' else color,markeredgewidth=1.3,
        linestyle='--' if name=='Re-FGO' else '-',
        label=name,zorder=4 if name=='Re-FGO' else 3)


def save(fig,out,name):
    for ext in ('pdf','svg','png'):
        fig.savefig(out/f'{name}.{ext}',dpi=240,bbox_inches='tight',pad_inches=.06)
    plt.close(fig)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    out=p.parse_args().output;out.mkdir(parents=True,exist_ok=True);appearance()
    data=generate_data(num_steps=100,emitter_radius=105,seed=7,gmm_weights=(.8,.2))
    results={}
    for mode in ('EKF','rEKF','iEKF','riEKF'):
        cfg=convert_kfv_to_fgo(Config(kfv=FilterConfig(mode=mode,max_iteration=10)))
        results[mode]=KfvEstimator(cfg,data).run()
        if mode=='riEKF':results['Re-FGO']=FgoEstimator(cfg,data).run()
    results['SW-FGO, W=5']=FgoEstimator(Config(fgo=FgoConfig(imitate_kfv=False,
        window_size=5,max_iteration=10,robust_kernel='huber')),data).run()
    order=list(STYLE);truth=data['true_positions'];time=np.arange(100)
    # Show the full sequence and a predeclared first 15 processed epochs.
    fig,axs=plt.subplots(1,2,figsize=(8.2,3.65),layout='constrained')
    for ax,ii,title in zip(axs,(np.arange(1,100),np.arange(1,16)),('(a) Complete trajectory','(b) First 15 updates')):
        ax.plot(*truth[:,ii],'-s',color='black',markersize=3,markevery=max(1,len(ii)//30),label='GT',zorder=2)
        for name in order:draw(ax,*results[name]['X_online'][:2,ii],name,every=max(1,len(ii)//30))
        axes(ax,'X (m)','Y (m)');ax.set_aspect('equal',adjustable='datalim');ax.set_title(title)
    # The actual extent on the right determines the outlined ROI.
    fig.canvas.draw();xmin,xmax=axs[1].get_xlim();ymin,ymax=axs[1].get_ylim()
    axs[0].add_patch(Rectangle((xmin,ymin),xmax-xmin,ymax-ymin,fill=False,color='k',linestyle='--',linewidth=1.2))
    handles,labels=axs[0].get_legend_handles_labels();fig.legend(handles,labels,loc='outside lower center',ncol=4,frameon=False)
    save(fig,out,'simulation_trajectory')
    # Current-epoch residual norm, not the graph's whole-window objective.
    fig,axs=plt.subplots(1,2,figsize=(8.2,3.3),layout='constrained')
    for name in order:
        x=results[name]['X_online'];err=np.linalg.norm(x[:2]-truth,axis=0)
        norm=np.linalg.norm(measurement_residuals(results[name],data),axis=0)
        draw(axs[0],time[1:],err[1:],name,3);draw(axs[1],time[1:],norm[1:],name,3)
    axes(axs[0],'Time (s)','2D position error (m)');axes(axs[1],'Time (s)','Range residual norm (m)')
    axs[0].set_title('(a) Online position error');axs[1].set_title('(b) Current-epoch residual')
    h,l=axs[0].get_legend_handles_labels();fig.legend(h,l,loc='outside lower center',ncol=3,frameon=False)
    save(fig,out,'simulation_error_residual')
    # A Gaussian nonlinear first update isolates the relinearization question.
    clean=generate_data(num_steps=100,emitter_radius=105,seed=7,gmm_weights=(1.,0.))
    cfg=convert_kfv_to_fgo(Config(kfv=FilterConfig(mode='iEKF',max_iteration=10)))
    filt,graph=KfvEstimator(cfg,clean).run(),FgoEstimator(cfg,clean).run()
    traces=[r['epoch_diagnostics'][1]['iteration_trace'] for r in (filt,graph)]
    fig,axs=plt.subplots(1,2,figsize=(8.2,3.6),layout='constrained')
    pts=np.array([t['state'][:2] for tr in traces for t in tr])
    xy=clean['true_positions'][:,1];lo=np.minimum(pts.min(0),xy);hi=np.maximum(pts.max(0),xy)
    padding=.12*np.maximum(1.,hi-lo);lo-=padding;hi+=padding
    xx,yy=np.meshgrid(np.linspace(lo[0],hi[0],150),np.linspace(lo[1],hi[1],150))
    rss=np.zeros_like(xx)
    for j,e in enumerate(clean['emitter_positions'].T):rss+=(clean['toa_measurements'][j,1]-np.hypot(xx-e[0],yy-e[1]))**2
    contour=axs[0].contour(xx,yy,np.sqrt(rss),levels=9,colors='#999999',linewidths=.7)
    axs[0].clabel(contour,inline=True,fontsize=10,fmt='%.0f')
    for tr,label,color,mark in zip(traces,('iEKF','Re-FGO'),(BLUE,RED),('x','o')):
        states=np.array([t['state'] for t in tr]);norm=[t['measurement_residual_norm_m'] for t in tr]
        kw=dict(color=color,marker=mark,markerfacecolor='none',label=label,markeredgewidth=1.4,markersize=5,linestyle='--' if label=='Re-FGO' else '-')
        axs[0].plot(states[:,0],states[:,1],**kw);axs[1].plot(np.arange(len(tr)),norm,**kw)
    axs[0].plot(*xy,'ks',markersize=5,label='GT')
    axes(axs[0],'X (m)','Y (m)');axes(axs[1],'Iteration','Range residual norm (m)')
    axs[0].set_title('(a) First-update path');axs[1].set_title('(b) Residual convergence')
    axs[0].legend(loc='lower left',framealpha=.95);axs[1].legend(frameon=False)
    save(fig,out,'simulation_iteration')
    rows=[]
    for mode in ('EKF','iEKF','rEKF','riEKF'):
        ep,ec=0.,0.
        for radius in (105.,1000.):
            for fraction in (0.,.2):
                d=generate_data(num_steps=100,emitter_radius=radius,seed=7,gmm_weights=(1-fraction,fraction))
                c=convert_kfv_to_fgo(Config(kfv=FilterConfig(mode=mode,max_iteration=10)))
                a,b=KfvEstimator(c,d).run(),FgoEstimator(c,d).run()
                ep=max(ep,float(np.max(np.abs(a['X_online'][:2,1:]-b['X_online'][:2,1:]))))
                ec=max(ec,float(np.max(np.abs(a['covariances'][1:,:2,:2]-b['covariances'][1:,:2,:2]))))
        rows.append(dict(mode=mode,max_position_difference_m=ep,max_position_covariance_difference_m2=ec))
    performance={name:dict(rmse_m=float(np.sqrt(np.mean(np.sum((r['X_online'][:2,1:]-truth[:,1:])**2,axis=0))))) for name,r in results.items()}
    metrics=dict(seed=7,epochs=100,dt_s=1,anchor_radius_m=105,initial_error_m=[100,-100],gmm_fraction=.2,
        gaussian_sigma_m=.1,outlier_sigma_m=10,iterations=10,huberdelta=2,comparison=rows,performance=performance)
    (out/'simulation_results.json').write_text(json.dumps(metrics,indent=2)+'\n')
    np.savez_compressed(out/'simulation_results.npz',truth=truth,measurements=data['toa_measurements'],
        **{name:r['X_online'] for name,r in results.items()})
    print(json.dumps(metrics,indent=2))

if __name__=='__main__':main()
