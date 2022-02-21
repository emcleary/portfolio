"""Example with streamed data

This example runs EKI by streaming a new model response into the
object on each iteration. This script demonstrates usage of the
initialize, run, and finalize methods. This approach is intended for
using with complex models that are computationally expensive and
cannot be written as simple functions.

This script should be run from the command line. Make sure `NumPy` is
installed in your Python environment.

"""


import sys
sys.path.append('../src')
sys.path.append('../tools')
import numpy as np
from eki import EKI
from models import NoisyLinear


if __name__=='__main__':

    n_iter = 20
    n_ens = 100
    n_var = 10
    n_par = 2

    # True parameters
    true_param = np.array([0.7, 2])

    # Generate parameter ensembles
    u_ens = np.random.uniform(0, 1, (true_param.size, n_ens))

    # Instantiate objects
    m = NoisyLinear(n_var, true_param)
    eki = EKI(u_ens, m.truth, m.cov)

    # EKI with streaming
    ui = eki.initialize_for_stream()
    for _ in range(n_iter):
        gi = m.run_model(ui)
        ui = eki.run_with_stream(gi)
    gi = m.run_model(ui)
    u_opt = eki.finalize_stream(gi)

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
