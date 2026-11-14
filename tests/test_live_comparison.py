from threading import Event
from contextlib import ExitStack
from unittest.mock import patch
import unittest
import numpy as np
from colab_comparison import prepare_comparison
from experiment_common import run_pair
from live_comparison import LiveComparison, run_interleaved
from core.objective import robust_weights, loss_cost


def case(kernel='huber', continuation=False, steps=30):
    return prepare_comparison('KFV vs FGO','Simulation',
        dict(method='riEKF',kernel=kernel,delta=1.,iterations=5,continuation=continuation),
        dict(method='Re-FGO'),simulation=dict(num_steps=steps),match_filter_settings=True)


class LiveComparisonTests(unittest.TestCase):
    def test_interleaving_preserves_batch_states_and_scrubbing(self):
        data,left,right=case()
        live=LiveComparison(data,[left,right]);events=[]
        def receive(i,e):
            events.append((e['epoch'],i));live.receive(i,e)
            if e['epoch']==10 and i==1:
                live.render();live.timeline.value=3
        streamed=run_interleaved(data,[left,right],receive,Event());live.finish()
        self.assertEqual(events,[(k,i) for k in range(1,30) for i in range(2)])
        self.assertEqual(live.timeline.value,3)
        self.assertFalse(live.follow.value)
        self.assertEqual(live.timeline.max,29)
        batch=run_pair(data,left,right)
        for a,b in zip(streamed['estimators'],batch['estimators']):
            np.testing.assert_allclose(a['result']['X_online'],b['result']['X_online'],atol=0,rtol=0)
        live.follow.value=True
        self.assertEqual(live.timeline.value,29)
        live.view.value='Residual norm'
        self.assertIn('Post-fit residual norm',live.svg)

    def test_stop_does_not_compute_remaining_epochs(self):
        data,left,right=case();stop=Event();events=[]
        def receive(i,event):
            events.append((i,event['epoch']))
            if i==1 and event['epoch']==4:stop.set()
        self.assertIsNone(run_interleaved(data,[left,right],receive,stop))
        self.assertEqual(len(events),8)

    def test_solver_time_excludes_waiting_for_live_view(self):
        data,left,right=case();clock=[0.]
        def now():
            clock[0]+=.001
            return clock[0]
        def receive(*_):clock[0]+=1000.
        with ExitStack() as stack:
            for module in ('Estimator','KfvEstimator','FgoEstimator'):
                stack.enter_context(patch(f'core.estimator.{module}.perf_counter',side_effect=now))
            experiment=run_interleaved(data,[left,right],receive,Event())
        for item in experiment['estimators']:
            self.assertGreater(item['result']['runtime_seconds'],0.)
            self.assertLess(item['result']['runtime_seconds'],1.)

    def test_redescending_continuation_preserves_matched_solvers(self):
        for kernel in ('cauchy','tukey'):
            with self.subTest(kernel=kernel):
                data,left,right=case(kernel,True,steps=100)
                experiment=run_interleaved(data,[left,right],lambda *_:None,Event())
                a,b=experiment['estimators']
                np.testing.assert_allclose(a['result']['X_online'],b['result']['X_online'],atol=1e-8,rtol=0)
                self.assertLess(a['metrics']['rmse'],2.)
                trace=a['result']['epoch_diagnostics'][1]['iteration_trace']
                self.assertGreater(trace[0]['robust_scale'],1.)
                self.assertEqual(trace[-1]['robust_scale'],1.)
                self.assertEqual(left['config'].kfv.robust_delta,1.)

    def test_tukey_has_exact_zero_weight_and_consistent_derivative(self):
        u=np.array([.2,.8,1.5,4.])
        for kernel in ('huber','cauchy','tukey'):
            for x,w in zip(u,robust_weights(u,kernel,1.)):
                derivative=(loss_cost(x+1e-6,kernel,1.)-loss_cost(x-1e-6,kernel,1.))/2e-6
                self.assertAlmostEqual(derivative,w*x,places=7)
        np.testing.assert_array_equal(robust_weights([2.,1e100],'tukey',1.),0.)


if __name__=='__main__':unittest.main()
