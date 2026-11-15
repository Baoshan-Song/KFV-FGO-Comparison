"""Three notebook workspaces with method-aware controls and shared data settings."""
import ipywidgets as W

from colab_comparison import MODULES, FILTERS, GRAPHS, prepare_comparison
from live_comparison import LiveComparison, run_interleaved
from threading import Thread, Event
from queue import SimpleQueue, Empty
import asyncio
from time import monotonic
from html import escape
from pathlib import Path
from datetime import datetime
from experiment_common import animate_payload, comparison_payload
from notebook_transport import ColabPanel, colab_output, encode_viewer, UI_REVISION


STYLE = '''<style>
.kfv-workspace {font-family:Arial,sans-serif;max-width:1160px;padding:16px;
background:#f4f7fb;border:1px solid #dce5ef;border-radius:14px;color:#223449}
.kfv-card {background:#fff;border:1px solid #dce5ef;border-radius:12px;padding:16px;margin:6px 0}
.kfv-card-blue {border-top:4px solid #1857d3}.kfv-card-red {border-top:4px solid #e3272b}
.kfv-heading {font-size:18px;font-weight:600;margin:0 0 5px}
.kfv-help {font-size:13px;line-height:1.5;color:#52677b;margin:0 0 9px}
.kfv-row {margin:5px 0;padding:5px 0;border-bottom:1px solid #eff3f7}
.kfv-icon {display:inline-block;background:#eaf2fb;color:#1857d3;border-radius:7px;
width:27px;text-align:center;padding:4px;margin-right:7px}
.kfv-label {font-size:14px;font-weight:600}.kfv-note {font-size:12px;color:#65788b;line-height:1.4}
.kfv-workspace .widget-label {color:#334b63}.kfv-workspace input:disabled,
.kfv-workspace select:disabled {color:#8592a0;background:#f0f3f6}
.kfv-workspace .widget-button {border-radius:7px}
</style>'''


def _heading(icon, title, description=''):
    return W.HTML(f'<div class="kfv-heading"><span class="kfv-icon">{icon}</span>{title}</div>'
                  f'<div class="kfv-help">{description}</div>')


def _row(icon, label, help_text, control):
    text = W.HTML(f'<span class="kfv-icon">{icon}</span><span class="kfv-label">{label}</span>'
                  f'<div class="kfv-note">{help_text}</div>', layout=W.Layout(width='44%', min_width='170px'))
    control.layout = W.Layout(width='54%', min_width='180px')
    return W.HBox([text, control], layout=W.Layout(align_items='center', flex_flow='row wrap',
                                                justify_content='space-between')).add_class('kfv-row')


def _card(children, width='100%'):
    return W.VBox(children, layout=W.Layout(width=width, min_width='300px')).add_class('kfv-card')


def _slider(value, low, high, step=1, integer=False):
    cls = W.IntSlider if integer else W.FloatSlider
    kwargs = {} if integer else {'readout_format':'.3g'}
    return cls(value=value, min=low, max=high, step=step, continuous_update=False, **kwargs)


def _log(value, low, high):
    return W.FloatLogSlider(value=value, base=10, min=low, max=high, step=.05,
                            readout_format='.3g', continuous_update=False)


class MethodCard:
    def __init__(self, kind, side, default):
        self.kind = kind
        self.method = W.Dropdown(options=FILTERS if kind == 'kfv' else GRAPHS, value=default)
        self.iterations = _slider(5, 1, 20, integer=True)
        self.kernel = W.Dropdown(options=['none','huber','cauchy','tukey'], value='huber')
        self.delta = _log(1., -2, 2)
        self._kernel_scales = {'none':1., 'huber':1., 'cauchy':2.385, 'tukey':4.685}
        self.window = _slider(5, 1, 50, integer=True)
        self.numerical = W.Checkbox(value=False, description='Numerical Jacobian', indent=False)
        self.continuation = W.Checkbox(value=True, description='Gradual scale', indent=False)
        self.continuation_row = _row('↘','Robust initialization','Starts broad, then decreases to the target over multiple updates. One-update rEKF uses the broad scale.',self.continuation)
        self.note = W.HTML()
        self._updating = False
        self._saved_iterations, self._saved_kernel, self._saved_window = 5, 'huber', 5
        self.iteration_row = _row('↻','Iterations','Maximum updates per epoch',self.iterations)
        self.kernel_row = _row('◇','Robust kernel','Reduce the influence of large whitened residuals',self.kernel)
        self.delta_row = _row('δ','Robust scale','Dimensionless scale, not a threshold in metres',self.delta)
        self.window_row = _row('▤','Window','States retained after a sliding-window solve',self.window)
        self.numerical_row = _row('∂','Derivative','Simulation only. Off uses the analytic Jacobian',self.numerical)
        rows = [_heading('A' if side == 'left' else 'B', f'{side.title()} method',
                         'Controls follow the selected estimator.'),
                _row('⇄','Estimator','Choose a filter variant' if kind == 'kfv' else 'Recursive or sliding-window graph',self.method),
                self.iteration_row,self.kernel_row,self.delta_row,self.continuation_row]
        if kind == 'fgo': rows += [self.window_row,self.numerical_row]
        rows += [self.note]
        self.widget = _card(rows,'49%').add_class('kfv-card-blue' if side == 'left' else 'kfv-card-red')
        self.method.observe(self._changed,names='value')
        self.kernel.observe(self._kernel_changed,names='value')
        self.apply_rules()

    def _changed(self, _):
        self.apply_rules()

    def _kernel_changed(self, change):
        if not self._updating:
            self._kernel_scales[change['old']] = self.delta.value
            self.delta.value = self._kernel_scales[change['new']]
            self.delta.disabled = self.kernel.disabled or self.kernel.value == 'none'
            self.continuation_row.layout.display = '' if self.kernel.value in ('cauchy','tukey') else 'none'

    def apply_rules(self, source='Simulation'):
        if self._updating:return
        self._updating=True
        try:
            if not self.iterations.disabled:self._saved_iterations=self.iterations.value
            if not self.kernel.disabled:self._saved_kernel=self.kernel.value
            if not self.window.disabled:self._saved_window=self.window.value
            single = self.kind == 'kfv' and self.method.value in ('EKF','rEKF')
            ordinary = self.kind == 'kfv' and self.method.value in ('EKF','iEKF')
            recursive = self.kind == 'fgo' and self.method.value == 'Re-FGO'
            self.iterations.value = 1 if single else self._saved_iterations
            self.iterations.disabled = single
            self.kernel.value = 'none' if ordinary else self._saved_kernel
            self.kernel.disabled = ordinary
            self.delta.disabled = ordinary or self.kernel.value == 'none'
            self.continuation_row.layout.display = '' if self.kernel.value in ('cauchy','tukey') else 'none'
            self.window.value = 1 if recursive else self._saved_window
            self.window.disabled = recursive
            self.numerical_row.layout.display = '' if source == 'Simulation' else 'none'
            text = ('EKF uses one update and no robust kernel.' if self.method.value == 'EKF' else
                    'IEKF enables repeated updates without robust weighting.' if self.method.value == 'iEKF' else
                    'rEKF uses one robust update.' if self.method.value == 'rEKF' else
                    'riEKF enables repeated robust updates.' if self.method.value == 'riEKF' else
                    'Re-FGO uses two QR marginalizations and retains one state.' if recursive else
                    'SW-FGO retains past states for relinearization.')
            self.note.value=f'<div class="kfv-note">{text}</div>'
        finally:self._updating=False

    def values(self):
        return dict(method=self.method.value,iterations=self.iterations.value,
                    kernel=self.kernel.value,delta=self.delta.value,window=self.window.value,
                    numerical_jacobian=self.numerical.value,continuation=self.continuation.value and self.kernel.value in ('cauchy','tukey'))


class ComparisonUI:
    def __init__(self, module):
        if module not in MODULES:raise ValueError(module)
        self.module=module
        self._busy=False
        self._colab_bridge=None
        self._previous={}
        self._closed=False
        self._viewer_revision=0
        self._viewer_packet=None
        defaults={'KFV vs KFV':('EKF','riEKF'),'KFV vs FGO':('iEKF','Re-FGO'),
                  'FGO vs FGO':('SW-FGO','SW-FGO')}
        self.left=MethodCard(MODULES[module][0],'left',defaults[module][0])
        self.right=MethodCard(MODULES[module][1],'right',defaults[module][1])
        if module=='FGO vs FGO':
            for card in (self.left,self.right):
                card.method.options=('SW-FGO',)
                card.method.disabled=True
            self.right.window.value=20
        self.source=W.ToggleButtons(options=['Simulation','Real GNSS/INS'],value='Simulation',
                                     icons=['flask','globe'],style={'button_width':'170px'})
        self.dataset=W.Dropdown(options=['medium','deep','harsh'],value='medium')
        self.dataset_row=_row('◉','Recording','Bundled UrbanNav data, loaded automatically',self.dataset)
        self.sim={
            'num_steps':_slider(100,30,300,10,True), 'anchor_radius':_slider(200,100,1000,25),
            'seed':W.IntText(value=7), 'white_sigma':_log(.1,-3,2),
            'outlier_weight':_slider(.2,0.,1.,.05), 'outlier_mean':_slider(0.,-100.,100.,1.),
            'outlier_sigma':_log(10.,-2,2),
            'initial_error_x':_slider(100.,-500.,500.,5.),
            'initial_error_y':_slider(-100.,-500.,500.,5.)}
        self.real={'dataset':self.dataset,'data_directory':W.Text(placeholder='Blank uses the selected recording'),
            'minimum_elevation_deg':_slider(-90.,-90.,90.,5.),'minimum_snr_dbhz':_slider(0.,0.,60.,1.),
            'position_process_std_m':_log(.3,-4,2),'velocity_process_std_mps':_log(.15,-4,2),
            'bias_process_std_mps2':_log(.01,-6,1),'clock_process_std_m':_log(100.,-2,4),
            'gps_minus_utc_s':W.BoundedIntText(value=18,min=0,max=100)}
        self.sim_data=_card([
            _row('●','Epochs','Number of states, including the initial prior',self.sim['num_steps']),
            _row('◎','Emitter radius (m)','Range-source geometry, not satellite altitude',self.sim['anchor_radius']),
            _row('σ','Range noise (m)','Nominal measurement standard deviation',self.sim['white_sigma'])])
        self.sim_faults=_card([
            _row('%','Outlier probability','Zero generates nominal Gaussian noise only',self.sim['outlier_weight']),
            _row('↥','Outlier mean (m)','Mean of the second Gaussian component',self.sim['outlier_mean']),
            _row('σ','Outlier spread (m)','Standard deviation of that component',self.sim['outlier_sigma'])])
        self.sim_initial=_card([
            _row('X','Initial X offset (m)','Shared by both methods',self.sim['initial_error_x']),
            _row('Y','Initial Y offset (m)','Shared by both methods',self.sim['initial_error_y']),
            _row('#','Random seed','Keep unchanged for a controlled comparison',self.sim['seed'])])
        self.sim_panel=W.Accordion(children=[self.sim_data,self.sim_faults,self.sim_initial],selected_index=None,
                                     titles=['◉ Geometry and nominal noise','⚑ Outlier simulation','↗ Initialization and repeatability'])
        self.real_noise=_card([
            _row('p','Position (m)','Process standard deviation per GNSS interval',self.real['position_process_std_m']),
            _row('v','Velocity (m/s)','Process standard deviation per GNSS interval',self.real['velocity_process_std_mps']),
            _row('b','Accelerometer bias (m/s²)','Process standard deviation per GNSS interval',self.real['bias_process_std_mps2']),
            _row('t','Receiver clock (m)','Process standard deviation per GNSS interval',self.real['clock_process_std_m'])])
        self.real_quality=_card([
            _row('∠','Minimum elevation (deg)','−90 disables the optional elevation gate',self.real['minimum_elevation_deg']),
            _row('▥','Minimum SNR (dB-Hz)','Zero disables the optional SNR gate',self.real['minimum_snr_dbhz'])])
        self.real_files=_card([
            _row('▣','Custom data directory','RINEX observations/navigation and Xsens data. Optional truth file.',self.real['data_directory']),
            _row('◷','GPS − UTC (s)','Historical offset for the recording',self.real['gps_minus_utc_s'])])
        self.real_panel=W.Accordion(children=[self.real_noise,self.real_quality,self.real_files],selected_index=None,
                                    titles=['σ Process uncertainty','◉ GNSS quality gates','▣ Custom data and time'])
        self.data_note=W.HTML()
        self.run=W.Button(description='Run comparison',icon='play',button_style='primary',layout=W.Layout(width='180px'))
        self.status=W.HTML('<span class="kfv-note">Ready. Plots update as each method computes.</span>')
        self.output=W.Output(layout=W.Layout(width='100%'))
        self.live_container=W.VBox()
        self.details=W.Button(description='Open full diagnostics',icon='line-chart',disabled=True,layout=W.Layout(width='210px'))
        self.details.layout.display='none'
        self.details.on_click(self._show_details)
        self.stop=W.Button(description='Stop',icon='stop',disabled=True,layout=W.Layout(width='90px'))
        self.stop_event=Event()
        self.use_cache=W.Checkbox(value=True,description='Use saved results',indent=False)
        self.stop.on_click(self._stop)
        self.connection=W.HTML(f'<div class="kfv-note">UI {UI_REVISION}</div>')
        self.widget=W.VBox([W.HTML(STYLE),_heading('⇄',module,'Choose data, configure the two methods, then run.'),
            _card([_heading('◉','Data source'),self.source,self.dataset_row,self.data_note]),
            W.HBox([self.left.widget,self.right.widget],layout=W.Layout(justify_content='space-between',flex_flow='row wrap')),
            _heading('⚙','Shared data settings','Optional. Defaults are ready to run.'),
            self.sim_panel,self.real_panel,
            W.HBox([self.run,self.stop,self.use_cache],layout=W.Layout(align_items='center',margin='14px 0 4px')),
            self.connection,self.status,
            W.HTML('<div class="kfv-note">Stop ends after the current epoch. Drag the live time slider to inspect completed results.</div>'),
            self.live_container,self.details,self.output],layout=W.Layout(width='100%')).add_class('kfv-workspace')
        self._syncing=False
        self.source.observe(self._sync,names='value')
        for card in (self.left,self.right):
            for control in (card.method,card.iterations,card.kernel,card.delta,card.continuation):control.observe(self._sync,names='value')
        self.run.on_click(self._run)
        self._sync()

    @property
    def match_filter_settings(self):
        return self.module=='KFV vs FGO' and self.right.method.value=='Re-FGO'

    def _sync(self, change=None):
        if self._syncing:return
        self._syncing=True
        try:
            real=self.source.value=='Real GNSS/INS'
            self.use_cache.layout.display='' if real else 'none'
            self.sim_panel.layout.display='none' if real else ''
            self.real_panel.layout.display='' if real else 'none'
            self.dataset_row.layout.display='' if real else 'none'
            self.data_note.value=('<div class="kfv-note">10-state pseudorange / IMU comparison with recorded AHRS attitude. Truth is used only for evaluation.</div>' if real else
                                  '<div class="kfv-note">Circular-motion range simulation. Both methods receive the same measurements and initial prior.</div>')
            self.left.apply_rules(self.source.value)
            self.right.apply_rules(self.source.value)
            if self.match_filter_settings:
                for key in ('iterations','kernel','delta','continuation'):
                    control=getattr(self.right,key);control.value=getattr(self.left,key).value;control.disabled=True
                self.right.numerical.disabled=True
                self.right.numerical.value=False
                self.right.note.value='<div class="kfv-note">Automatically follows the KFV settings on the left. Re-FGO is solved independently.</div>'
            else:
                self.right.numerical.disabled=False
                self.right.continuation.disabled=False
        finally:self._syncing=False

    def _show_details(self, _):
        if not hasattr(self,'experiment') and not hasattr(self,'cached_payload'):return
        payload=(self.cached_payload if hasattr(self,'cached_payload')
                 else comparison_payload(self.experiment))
        payload=dict(payload,compact=True)
        self._viewer_packet=encode_viewer(payload)
        self._viewer_revision+=1
        if self._live is not None:self._live.play.playing=False
        # Keep the existing result visible until the browser confirms the old
        # renderer has actually produced an SVG. Sending data is not rendering.
        self.output.outputs=()
        if self._colab_bridge is None:
            self.output.append_display_data(animate_payload(payload,interval=80))

    def _viewer_ready(self, revision):
        if revision==self._viewer_revision and self._viewer_packet is not None:
            self.live_container.layout.display='none'

    def _stop(self, _=None):
        if not self._busy:return
        self.stop_event.set()
        self._stage='Stopping after the current operation…'
        self.status.value=f'<b>{self._stage}</b>'

    def _run(self, _):
        if self._busy or self.run.disabled:return
        try:
            # Colab explicitly polls from fresh kernel callbacks. Do not depend
            # on an asyncio task surviving the widget's originating request.
            loop=None if self._colab_bridge else asyncio.get_running_loop()
            self._busy=True
            self._finished=False
            if getattr(self,'_live',None) is not None:self._live.play.playing=False
            self._live=None
            self._started=monotonic()
            self._stage='Starting worker…'
            self._events=SimpleQueue()
            self.run.disabled=True;self.stop.disabled=False;self.details.disabled=True
            self.stop_event.clear()
            self.status.value=f'<b>{self._stage}</b>'
            # Direct state updates avoid output capture / clear_output messages
            # tied to a completed cell or widget callback.
            self.output.outputs=()
            self._viewer_packet=None
            self._viewer_revision+=1
            self.live_container.children=[]
            self.live_container.layout.display=''
            for field in ('experiment','cached_payload'):
                if hasattr(self,field):delattr(self,field)
            args=(self.module,self.source.value,self.left.values(),self.right.values())
            options=dict(simulation={k:v.value for k,v in self.sim.items()},
                         real={k:v.value for k,v in self.real.items()},match_filter_settings=self.match_filter_settings)
            reuse=self.use_cache.value
            controls=[self.source,self.dataset,self.use_cache,*self.sim.values(),*self.real.values()]
            for card in (self.left,self.right):
                controls += [card.method,card.iterations,card.kernel,card.delta,card.window,card.numerical,card.continuation]
            self._previous={id(c):(c,c.disabled) for c in controls}
            for control,_ in self._previous.values():control.disabled=True
            self._start_worker(args,options,reuse)
            if loop is not None:self.ui_task=loop.create_task(self._pump())
        except Exception as error:
            self.stop_event.set()
            self._stage=f'Start failed: {type(error).__name__}: {error}'
            self.status.value=f'<b>{escape(self._stage)}</b>'
            self._finished=True
            if not getattr(self,'worker',None) or not self.worker.is_alive():self._restore_controls()

    def _start_worker(self,args,options,reuse):
        def report(message):
            if self.stop_event.is_set():raise InterruptedError('Stopped')
            self._events.put(('status',message))
        def worker():
            try:
                if args[1]=='Real GNSS/INS' and reuse:
                    from run_cache import load_result_cache
                    report('Checking saved results for these exact settings')
                    _,left,right=prepare_comparison(*args,**options,load_observations=False)
                    payload=load_result_cache([left,right])
                    report('Saved-result check complete')
                    if payload is not None:
                        self._events.put(('cached',payload));return
                data,left,right=prepare_comparison(*args,**options,progress=report,stop=self.stop_event)
                if self.stop_event.is_set():
                    raise InterruptedError('Stopped during loading')
                self._events.put(('ready',(data,[left,right])))
                experiment=run_interleaved(data,[left,right],lambda i,e:self._events.put(('epoch',(i,e))),self.stop_event)
                if experiment is None:
                    raise InterruptedError('Stopped. Completed epochs remain available.')
                report('Saving completed results')
                if args[1]=='Real GNSS/INS':
                    from real_experiment import save_real_experiment
                    from run_cache import save_result_cache
                    folder=Path('artifacts/colab')/self.module.lower().replace(' ','-')/datetime.now().strftime('%Y%m%d-%H%M%S-%f')
                    save_real_experiment(experiment,folder)
                    save_result_cache(experiment)
                self._events.put(('complete',experiment))
            except InterruptedError as error:
                self._events.put(('stopped',str(error)))
            except Exception as error:
                self._events.put(('error',f'{type(error).__name__}: {error}'))
            finally:self._events.put(('done',None))
        self.worker=Thread(target=worker,name='kfv-comparison',daemon=True)
        self.worker.start()

    def _drain_events(self):
        """A bounded, synchronous UI update, also callable without asyncio."""
        if not self._busy:return
        try:
            deadline=monotonic()+.035
            while not self._finished and monotonic()<deadline:
                try:kind,value=self._events.get_nowait()
                except Empty:break
                if kind=='status':
                    if not self.stop_event.is_set():self._stage=value
                elif kind=='ready':
                    self._live=LiveComparison(*value)
                    self.live_container.children=[self._live.widget]
                    if not self.stop_event.is_set():self._stage='Computing both methods'
                elif kind=='epoch':self._live.receive(*value)
                elif kind=='cached':
                    if self.stop_event.is_set():
                        self._stage='Stopped'
                        continue
                    self.cached_payload=value
                    self._live=LiveComparison.from_payload(value)
                    self.live_container.children=[self._live.widget]
                    self.details.disabled=False
                    self._show_details(None)
                    self._stage=f'Loaded saved results · {value["frames"]-1} updates · stored solver timings'
                elif kind=='complete':
                    self.experiment=value;self._live.finish();self.details.disabled=False
                    self._show_details(None)
                    self._stage='Complete. Play or scrub the timeline below'
                elif kind in ('stopped','error'):
                    if self._live is not None:self._live.finish(stopped=True)
                    self._stage=value
                    if kind=='error':self.output.append_stdout(value+'\n')
                elif kind=='done':self._finished=True
            self.status.value=f'<b>{escape(self._stage)}</b> · {monotonic()-self._started:.1f} s'
        except Exception as error:
            self.stop_event.set()
            self._finished=True
            self._stage=f'Display error: {type(error).__name__}: {error}'
            self.status.value=f'<b>{escape(self._stage)}</b>'
        if self._finished and not self.worker.is_alive():self._restore_controls()

    def _restore_controls(self):
        for control,disabled in self._previous.values():control.disabled=disabled
        self._previous={}
        self._busy=False
        self.run.disabled=False;self.stop.disabled=True;self._sync()

    async def _pump(self):
        # Non-Colab notebooks use their event loop. Colab takes the explicit
        # request path instead (see notebook_transport.py).
        while self._busy:
            self._drain_events()
            if self._busy:await asyncio.sleep(.05)


def comparison_ui(module):
    """Return the workspace, with an explicit request transport in Colab."""
    ui=ComparisonUI(module)
    output=colab_output()
    return ColabPanel(ui,output) if output is not None else ui.widget
