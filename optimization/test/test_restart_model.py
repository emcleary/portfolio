"""Test case: restarting the model implementation

This test ensures that an EKI object using the model can be restarted
properly.

"""


import pytest
import sys
sys.path.append('../src')
sys.path.append('../tools')
import numpy as np
from eki import EKI
from models import Model


def test_main():

    # Seed for RNG
    if __name__ == '__main__':
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
    eki_1 = EKI(u_ens, m.truth, m.cov, model=m.run_model)
    eki_2 = EKI(u_ens, m.truth, m.cov, model=m.run_model)

    # EKI with the model
    np.random.seed(seed)
    u_opt_1 = eki_1.run_with_model(2 * n_iter)

    # EKI with streaming
    np.random.seed(seed)
    _ = eki_2.run_with_model(n_iter)
    u_opt_2 = eki_2.run_with_model(n_iter)

    # Check errors
    assert np.all(eki_1.error == eki_2.error), "Errors are not identical"

    # Check optimized parameters
    assert np.all(u_opt_1 == u_opt_2), "Optimized parameters are not identical"


if __name__=='__main__':
    test_main()
