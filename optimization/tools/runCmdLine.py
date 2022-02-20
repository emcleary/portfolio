import sys
sys.path.append('../src')
import numpy as np
from eki import EKI
from argparse import ArgumentParser
from models import MODEL

def get_cmd_line_arguments():
    parser = ArgumentParser(description="Experiment with optimizing a noisy, linear model")

    parser.add_argument('n_iter', type=int, help='Number of iterations for optimizing')
    parser.add_argument('-n_ens', type=int, default=100, help='Ensemble size')
    parser.add_argument('-n_var', type=int, default=5, help='Dimensionality of the data space')
    parser.add_argument('-n_par', type=int, default=2, help='Dimensionality of the parameter space')
    parser.add_argument('-min_param', type=int, default=0, help='Minimum value of the true parameters')
    parser.add_argument('-max_param', type=int, default=1, help='Maximum value of the true parameters')
    parser.add_argument('-seed', type=int, default=None, help='Seed for the random number generator')

    return parser.parse_args()


if __name__=='__main__':

    args = get_cmd_line_arguments()
    n_iter = args.n_iter
    n_ens = args.n_ens
    n_var = args.n_var
    n_par = args.n_par
    min_param = args.min_param
    max_param = args.max_param
    seed = args.seed

    if n_par > n_var:
        print('The number of model parameters cannot be larger than the dimensionality of the data.\n'
              'Otherwise there is not a unique solution to the optimized parameters.')

    if seed is not None:
        np.random.seed(seed)

    # True parameters
    true_param = np.random.uniform(min_param, max_param, n_par)

    # Generate parameter ensembles
    u_ens = np.random.uniform(min_param, max_param, (n_par, n_ens))

    # Optimize model parameters
    m = MODEL(n_var, true_param)
    eki = EKI(u_ens, m.truth, m.cov, model=m.run_model)
    u_opt = eki.run_with_model(u_ens, n_iter)

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
