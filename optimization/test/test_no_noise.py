import pytest
import sys
sys.path.append('../src')
sys.path.append('../tools')
import numpy as np
from eki import EKI
from models import Model


def test_main():

    np.random.seed(42)

    # Ensemble size
    n_ens = 100

    # Truth
    n_var = 5
    true_param = np.array([0.7, 2])

    # Generate ensembles
    u_ens = np.random.normal(loc=0, scale=2, size=(true_param.size, n_ens))

    # Model and EKI objects
    m = Model(n_var, true_param)
    eki = EKI(u_ens, m.truth, m.cov, model=m.run_model)

    # EKI iteration
    u_opt = eki.run_with_model(1)

    # Converges in one iterations
    err = np.abs(true_param - eki.u[-1].mean(1))
    assert all([ei < 1e-3 for ei in err])


if __name__=='__main__':
    test_main()
