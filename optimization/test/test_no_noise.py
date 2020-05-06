import pytest
import sys
sys.path.append('../src')
import numpy as np
from eki import EKI

class MODEL:

    # numvar: number of variables
    def __init__(self, numvar, parameters):
        # True parameters
        self.p = parameters
        numpar = parameters.size

        # Linear model
        self.A = np.random.normal(loc=0, scale=2, size=(numvar, numpar))

        # Covariance for noise
        self.cov = np.zeros((numvar,numvar))

        # True model
        self.truth = self.f(self.p)

    # Model without noise
    def f(self, parameters):
        return self.A.dot(parameters)


def test_main():

    np.random.seed(42)
    
    # Ensemble size
    J = 100

    # Truth
    numvar = 5
    parameters = np.array([0.7, 2])
    m = MODEL(numvar, parameters)

    # Generate ensembles
    u_ens = np.random.normal(loc=0, scale=2, size=(parameters.size, J))

    # Ensemble Kalman Inversion object
    eki = EKI(u_ens, m.truth, m.cov)
    g_ens = m.f(u_ens)

    # EKI iteration
    eki.update(g_ens)

    # Converges in one iterations
    err = np.abs(parameters - eki.u[-1].mean(1))
    assert all([ei < 1e-14 for ei in err])
