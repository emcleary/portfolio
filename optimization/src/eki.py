import numpy as np
import matplotlib.pyplot as plt

class EKI:

    # y = g(u) + N(0, cov)
    def __init__(self, parameters, truth, cov):

        # Truth statistics
        self.g_t = truth
        self.cov = cov

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

        return


    # Compute error using mean of ensemble model evaluations.
    def compute_error(self, iteration=-1):
        diff = self.g_t - self.g[iteration].mean(1)
        error = diff.dot(np.linalg.solve(self.cov, diff))
        self.error = np.append(self.error, error)


    # g: data, i.e. g(u), with shape (num_ensembles, num_elements)
    def update(self, g):

        u = np.copy(self.u[-1])

        # Ensemble statistics
        u_bar = u.mean(1)
        p_bar = g.mean(1)
        c_up = np.tensordot(u, g.T, 1) / self.J
        c_pp = np.tensordot(g, g.T, 1) / self.J

        # Update parameter ensemble
        y = self.g_t + np.random.multivariate_normal(np.zeros(self.g_t.size), self.cov, size=self.J) 
        tmp = np.linalg.solve(c_pp + self.cov, y.T-g)
        u += c_up.dot(tmp)

        # Store parameters and observations
        self.u = np.append(self.u, [u], axis=0)
        self.g = np.append(self.g, [g], axis=0)

        return
