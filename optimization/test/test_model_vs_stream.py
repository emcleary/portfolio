"""Test case: model and streaming implementations

This test ensures that both model and streaming methods in the EKI
class will return identical results.

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
    eki_model = EKI(u_ens, m.truth, m.cov, model=m.run_model)
    eki_stream = EKI(u_ens, m.truth, m.cov)

    # EKI with the model
    np.random.seed(seed)
    u_opt_model = eki_model.run_with_model(n_iter)

    # EKI with streaming
    np.random.seed(seed)
    ui = eki_stream.initialize_for_stream()
    for _ in range(n_iter):
        gi = m.run_model(ui)
        ui = eki_stream.run_with_stream(gi)
    gi = m.run_model(ui)
    u_opt_stream = eki_stream.finalize_stream(gi)

    # Check errors
    assert np.all(eki_model.error == eki_stream.error), "Errors are not identical"

    # Check optimized parameters
    assert np.all(u_opt_model == u_opt_stream), "Optimized parameters are not identical"


if __name__=='__main__':
    test_main()
