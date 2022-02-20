import numpy as np
from typing import Optional, Generic
from myTypes import *

class EKI(Generic[NDArrayFloat]):
    p: int
    J: int
    g_t: NDArrayFloat
    cov: NDArrayFloat
    cov_inv: NDArrayFloat
    u: NDArrayFloat
    error: NDArrayFloat
    run_Model: Optional[CallModel]

    # y = g(u) + N(0, cov)
    def __init__(self, parameters: NDArrayFloat, truth: NDArrayFloat,
                 cov: NDArrayFloat, model: Optional[CallModel] = None) -> None:

        # Parameter (p) and ensemble sizes (J)
        assert parameters.ndim == 2, 'parameters must be a 2-dimensional array (num parameters, ensemble)'
        self.p, self.J = parameters.shape

        # Check other inputs
        assert truth.ndim == 1, 'truth must be a 1-dimensional array'
        assert cov.ndim == 2, 'cov must be a 2-dimensional array'

        # Truth statistics
        self.g_t = truth
        self.cov = cov

        # Ensure invertibility of the covariance
        # (typically only a problem here on models with no noise)
        try:
            self.cov_inv = np.linalg.inv(self.cov)
        except np.linalg.LinAlgError:
            self.cov += 1e-8 * np.eye(cov.shape[0], dtype=self.cov.dtype)
            self.cov_inv = np.linalg.inv(self.cov)

        # Parameters with shape (iteration, number of parameters, ensemble size)
        self.u = parameters[np.newaxis]

        # Store observations during iterations
        self.g = np.empty((0,truth.size,self.J))

        # Error
        self.error = np.empty(0)

        # Model to be trained
        self.run_model = model


    # Compute error using mean of ensemble model evaluations.
    def _compute_error(self, iteration: int = -1):
        assert len(self.u) == len(self.g), 'Iterations do not match up between the model parameters and model respone.'
        diff = self.g_t - self.g[iteration].mean(1)
        error = diff @ self.cov_inv @ diff
        self.error = np.append(self.error, error)


    def update_parameters(self) -> None:

        u = np.copy(self.u[-1])
        g = np.copy(self.g[-1])

        # Ensemble statistics
        c_up = np.tensordot(u, g.T, 1) / self.J
        c_pp = np.tensordot(g, g.T, 1) / self.J

        # Update parameter ensemble
        y = self.g_t + np.random.multivariate_normal(np.zeros(self.g_t.size), self.cov, size=self.J) 
        tmp = np.linalg.solve(c_pp + self.cov, y.T-g)
        u += c_up @ tmp

        # Store parameters and observations
        self.u = np.append(self.u, [u], axis=0)


    # g is the evaluation of a model with shape (num_ensembles, num_elements)
    def update_model_response(self, g: NDArrayFloat) -> None:
        assert g.shape == (self.g_t.size, self.J), "g must have the shape of data dimensionality x ensemble size"
        self.g = np.append(self.g, [g], axis=0)
        self._compute_error()


    def run_with_model(self, n_iterations: int) -> NDArrayFloat:
        assert self.run_model, "EKI has not been given a model to use!"
        g_ens = self.run_model(self.u[-1])
        self.update_model_response(g_ens)
        for _ in range(n_iterations):
            self.update_parameters()
            g_ens = self.run_model(self.u[-1])
            self.update_model_response(g_ens)
        return self.u[-1].mean(1)
