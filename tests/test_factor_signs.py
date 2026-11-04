import unittest
import numpy as np
from core.fgo.factor import State, PdrFactor, AutoDiffFactor, AutoDiffRobustFactor

class FactorSignTests(unittest.TestCase):
    def test_finite_difference_residuals_point_towards_zero(self):
        for factor in (
            AutoDiffFactor([State(1,1,np.zeros(4))], np.ones(4), np.eye(4), lambda x,z:z-x),
            AutoDiffRobustFactor([State(1,1,np.zeros(4))], np.ones(4), lambda x,z:z-x, lambda e:e)):
            factor.evaluate()
            delta=np.linalg.lstsq(factor.A,factor.b,rcond=None)[0]
            np.testing.assert_allclose(delta, np.ones(4), atol=1e-10)

    def test_displacement_factor_uses_second_minus_first(self):
        old,new=State(1,1,np.zeros(4)),State(2,2,np.zeros(4))
        f=PdrFactor([old,new],np.ones(4),np.eye(4)).evaluate()
        delta=np.linalg.solve(f.A[:,4:],f.b)
        new.value+=delta
        np.testing.assert_allclose(f.evaluate().b,0.)
