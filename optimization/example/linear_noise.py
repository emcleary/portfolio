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
        C = np.random.normal(loc=0, scale=2, size=(numvar, numvar))
        self.cov = 0.2*(C.T.dot(C) + 0.1*np.eye(numvar))

        # True model
        self.truth = self._ft(self.p)

    # Model without noise
    def _ft(self, parameters):
        return self.A.dot(parameters)

    # Noisy model evaluation
    # parameters: (num variables, ensemble size)
    def f(self, parameters):
        if parameters.ndim == 1:
            parameters = np.expand_dims(parameters, axis=1)

        # Evaluate model
        g = self._ft(parameters)

        # Evaluate noise
        numvar, numens = g.shape
        eta = np.random.multivariate_normal(np.zeros(numvar), self.cov, size=numens).T

        return g + eta


def main():

    # Maximum iterations
    iter_max = 1000

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

    # Iterate
    iter = 0
    while iter < iter_max:
        eki.update(g_ens)
        eki.compute_error()
        if eki.error[-1] < 0.01: break
        iter += 1
        u_ens = eki.u[-1]
        g_ens = m.f(u_ens)

    print('')
    if iter == 1:
        print(iter,'iteration')
    else:
        print(iter,'iterations')

    print('')
    print('Mean parameters')
    print(u_ens.mean(1))
    print('')
    print('True parameters')
    print(m.p)
    print('')
    print('Error')
    print(eki.error[-1])
    print('')


main()
