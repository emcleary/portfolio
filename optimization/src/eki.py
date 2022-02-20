import numpy as np
import matplotlib.pyplot as plt

class EKI:

    # y = g(u) + N(0, cov)
    def __init__(self, parameters, truth, cov, model=None):

        # Truth statistics
        self.g_t = truth
        self.cov = cov

        # Ensure invertibility of the covariance
        # (typically only a problem here on models with no noise)
        try:
            self.cov_inv = np.linalg.inv(self.cov)
        except np.linalg.LinAlgError:
            self.cov += 1e-8 * np.eye(cov.shape[0])
            self.cov_inv = np.linalg.inv(self.cov)

        # Parameters
        self.u = parameters[np.newaxis]

        # Ensemble size
        self.p, self.J = parameters.shape

        # Size of statistics
        self.n_obs = truth.size

        # Store observations during iterations
        self.g = np.empty((0,truth.size,self.J))

        # Error
        self.error = np.empty(0)

        # Model to be trained
        self.run_model = model


    # Compute error using mean of ensemble model evaluations.
    def compute_error(self, iteration=-1):
        diff = self.g_t - self.g[iteration].mean(1)
        error = diff @ self.cov_inv @ diff
        self.error = np.append(self.error, error)


    def update_parameters(self):

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
    def update_model_response(self, g):
        assert g.shape == (self.g_t.size, self.J), "g must have the shape of data dimensionality x ensemble size"
        self.g = np.append(self.g, [g], axis=0)


    def run_with_model(self, param_ens, n_iterations):
        assert self.run_model, "EKI has not been given a model to use!"
        g_ens = self.run_model(param_ens)
        self.update_model_response(g_ens)
        for _ in range(n_iterations):
            self.compute_error()
            self.update_parameters()
            g_ens = self.run_model(self.u[-1])
            self.update_model_response(g_ens)
        self.compute_error()
        return self.u[-1].mean(1)
