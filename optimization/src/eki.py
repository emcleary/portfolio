import numpy as np
from typing import Optional, Generic
from myTypes import *

class EKI(Generic[NDArrayFloat]):
    """
    Optimizes model parameters with Ensemble Kalman Inversion (EKI).

    Attributes
    ----------
    p : int
        The number of model parameters
    J : int
        The ensemble size
    g_t : 1d array
        True data uses for training model parameters
    cov : 2d array
        Covariance of the (true) data
    cov_inv : 2d array
        Inverse of the covariance used for error norms
    u : 3d array
        Stores ensembles of parameters for all iterations
    g : 3d array
        Stores ensembles of model responses for all iterations
    error : 1d array
        Errors calculated from model responses and the truth
    run_model : function, optional
        A function for evaluating the model

    Methods
    -------
    get_parameters()
        Returns and ensemble average of the last iteration of parameters.
    run_with_model(N)
        Runs the EKI algorithm using a model function for N iterations. 
    initialize_for_stream()
        Ensures the object is ready for streaming iterations and 
        returns model parameter ensemble.
    run_with_stream(g)
        Runs 1 iteration of EKI using model response g, returning the 
        updated model parameter ensemble.
    finalize_stream(g)
        Finalizes error calculations with g and returns the optimized 
        ensemble average of the parameters.
    """
    
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
        """
        Parameters
        ----------
        parameters : 2d array
            An ensemble of model parameters, used for defining the 
            number of parameters and the ensemble size, as well as
            for the first iteration of EKI.
        truth : 1d array
            True, typically observational data used to train the model.
        cov : 2d array
            A covariance of the data space, typically calculated from
            observational data or estimated.
        model : function, optional
            A function taking in an ensemble of model parameters and
            used to generate data, if needed.
        """

        # Parameter (p) and ensemble sizes (J)
        assert parameters.ndim == 2, 'parameters must be a 2-dimensional array (num parameters, ensemble)'
        self.p, self.J = parameters.shape

        # Check other inputs
        assert truth.ndim == 1, 'truth must be a 1-dimensional array'
        assert cov.ndim == 2, 'cov must be a 2-dimensional array'
        assert cov.shape == (truth.size, truth.size), 'covariance and truth must represent the same number of variables'

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


    def get_parameters(self) -> NDArrayFloat:
        """
        Calculates and returns and ensemble average of parameters.

        Returns
        -------
        array
            A 1d array representing the ensemble average of model
            parameters from the most recent iteration.
        """
        return self.u[-1].mean(1)


    def _update_parameters(self) -> None:
        """
        The EKI algorithm.

        Reference: 
            Marco A Iglesias et al., "Ensemble Kalman methods for 
            inverse problems". 2013, Inverse Problems, 29, 045001
        """
        assert len(self.u) == len(self.g), "Iterations for parameters and model responses do not match!"

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


    def _update_model_response(self, g: NDArrayFloat) -> None:
        assert g.shape == (self.g_t.size, self.J), "g must have the shape of data dimensionality x ensemble size"
        assert len(self.u) == len(self.g) + 1, "Must have 1 more iteration of model parameters than model responses."
        self.g = np.append(self.g, [g], axis=0)
        self._compute_error()


    def run_with_model(self, n_iterations: int) -> NDArrayFloat:
        """
        Run EKI, modeling data on the fly using the run_model function.

        Parameters
        ----------
        n_iterations : int
            Number of iterations of the EKI algorithm

        Returns
        -------
        array
            An ensemble average of model parameters from the last
            iteration.
        """
        
        assert self.run_model, "EKI has not been given a model to use!"

        # Evaluate the model as needed.
        # Should only matter the first time this function is called.
        if len(self.u) == len(self.g) + 1:
            g_ens = self.run_model(self.u[-1])
            self._update_model_response(g_ens)

        # EKI algorithm
        for _ in range(n_iterations):
            self._update_parameters()
            g_ens = self.run_model(self.u[-1])
            self._update_model_response(g_ens)

        # Return optimized parameters
        return self.get_parameters()


    def initialize_for_stream(self) -> NDArrayFloat:
        """
        Ensures the class instance is setup for EKI iterations with
        streamed data.

        Returns
        -------
        array
            A 2d array of the most recent ensemble of model parameters.
        """
        if len(self.u) == len(self.g):
            self._update_parameters()
        return self.u[-1]


    def run_with_stream(self, g: NDArrayFloat) -> NDArrayFloat:
        """
        Run one iteration of EKI with the supplied ensemble of model
        responses.

        Parameters
        ----------
        g : array
            A 2d array containing an ensemble of model responses
            calculated using the most recent ensemble of parameters.

        Returns
        -------
        array
            A 2d array of the updated ensemble of model parameters.
        """
        self._update_model_response(g)
        self._update_parameters()
        return self.u[-1]


    def finalize_stream(self, g: NDArrayFloat) -> NDArrayFloat:
        """
        Add the final model response and return tuned parameters.

        Returns
        -------
        array
            A 1d array of optimized model parameters.
        """
        self._update_model_response(g)
        return self.get_parameters()
