"""Example with a linear, noisy model

This example runs EKI by passing the model into the class.

This script should be run from the command line. Make sure `NumPy` is
installed in your Python environment.

"""


import sys
sys.path.append('../src')
sys.path.append('../tools')
import numpy as np
from eki import EKI
from models import NoisyLinear


def main():

    n_iter = 20
    n_ens = 100
    n_var = 10
    n_par = 2

    # True parameters
    true_param = np.array([0.7, 2])

    # Generate parameter ensembles
    u_ens = np.random.uniform(0, 1, (true_param.size, n_ens))

    # Optimize model parameters
    m = NoisyLinear(n_var, true_param)
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


if __name__=='__main__':
    main()
