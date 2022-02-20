import sys
sys.path.append('../src')
sys.path.append('../tools')
import numpy as np
from eki import EKI
from models import Model
from myTypes import NDArrayFloat64

class Nonlinear(Model):
    def initialize_model(self) -> None:
        self.A = np.random.normal(0, 2, size=(self.n_variables, self.true_param.size))
        self.B = np.random.normal(0, 2, size=(self.n_variables, self.true_param.size))
        C = np.random.normal(loc=0, scale=2, size=(self.n_variables, self.n_variables))
        self.cov = 0.01*(C.T.dot(C) + 0.1*np.eye(self.n_variables))


    def run_ensemble(self, parameters: NDArrayFloat64) -> NDArrayFloat64:
        g = self.A @ parameters + 0.5 * (self.B @ parameters) ** 2
        n_variables, n_particles = g.shape
        eta = np.random.multivariate_normal(np.zeros(n_variables), self.cov, size=n_particles).T
        return g + eta
    


if __name__=='__main__':

    n_iter = 20
    n_ens = 100
    n_var = 5
    n_par = 1

    # True parameters
    true_param = np.random.uniform(0, 1, n_par)

    # Generate parameter ensembles
    u_ens = np.random.uniform(0, 1, (true_param.size, n_ens))

    # Optimize model parameters
    m = Nonlinear(n_var, true_param)
    eki = EKI(u_ens, m.truth, m.cov, model=m.run_model)
    u_opt = eki.run_with_model(n_iter)

    print('')
    print('Optimized parameters')
    print(u_opt)
    print('')
    print('True parameters')
    print(true_param)
    print('')
    print('Error')
    print(eki.error[-1])
    print('')
