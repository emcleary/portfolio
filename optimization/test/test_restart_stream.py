"""Test case: restarting the streaming implementation

This test ensures that an EKI object using the streaming methods can
be restared properly.

"""


import pytest
import sys
sys.path.append('../src')
sys.path.append('../tools')
import numpy as np
from eki import EKI
from models import Model


def streaming(eki, model, n_iter):

    # Make sure parameters and model responses are up to date.
    # NB: this does 1 iteration of EKI when restarting.
    ui = eki.initialize_for_stream()

    # Go back and forth between running the EKI algorithm and the
    # model.
    for _ in range(n_iter):
        gi = model.run_model(ui)
        ui = eki.run_with_stream(gi)
    gi = model.run_model(ui)

    # Finalize the error calculations and return the optimized model
    # parameters.
    return eki.finalize_stream(gi)


def test_main(random_seed=False):

    # Seed for RNG
    if random_seed:
        seed = np.random.randint(1, 2**31-1)
        print('Seed set to', seed)
    else:
        seed = 42

    # Number of EKI iterations
    n_iter = 10
    
    # Ensemble size
    n_ens = 100

    # Truth
    n_var = 5
    true_param = np.array([0.7, 2])

    # Generate ensembles
    u_ens = np.random.normal(loc=0, scale=2, size=(true_param.size, n_ens))

    # Model and streaming objects
    m = Model(n_var, true_param)
    eki_1 = EKI(u_ens, m.truth, m.cov)
    eki_2 = EKI(u_ens, m.truth, m.cov)

    # EKI with the model
    np.random.seed(seed)
    u_opt_1 = streaming(eki_1, m, 2 * n_iter)

    # EKI with streaming
    np.random.seed(seed)
    _ = streaming(eki_2, m, n_iter)
    u_opt_2 = streaming(eki_2, m, n_iter-1)

    # Check errors
    assert np.all(eki_1.error == eki_2.error), "Errors are not identical"

    # Check optimized parameters
    assert np.all(u_opt_1 == u_opt_2), "Optimized parameters are not identical"


if __name__=='__main__':
    test_main(True)
