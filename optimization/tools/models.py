import sys
sys.path.append('../src')
import numpy as np
from argparse import ArgumentParser
from typing import Optional, Generic
from my_types import *


class Model(Generic[NDArrayFloat]):
    """
    A model template

    Attributes
    ----------
    true_param : array
        A 1d array of target parameters set by the user.
    n_dim : int
        The dimenion of the data space.
    truth : array
        The model response evaluated at the true parameters.
    cov : array
        The covariance of the model response.

    Methods
    -------
    initialize_model()
        Generate or set all hyperparameters for the model.
    run_ensemble(param_ens)
        Run the model on an ensemble of parameters.
    run_model(param)
        Run the model on either a set of parameters or ensemble of
        parameters.
    """
    
    true_param: NDArrayFloat
    n_dim: int
    _truth: Optional[NDArrayFloat]
    _cov: Optional[NDArrayFloat]

    # n_dim: number of variables
    # parameters: numpy array of true parameters
    def __init__(self, n_dim: int, true_param: NDArrayFloat) -> None:
        """
        Parameters
        ----------
        n_dim : int
            The dimension of the data space
        true_param : array
            A 1d array of target parameters.
        """
        self.true_param = true_param
        self.n_dim = n_dim

        self._truth = None
        self._cov = None
        self.initialize_model()


    @property
    def truth(self) -> NDArrayFloat:
        if self._truth is None:
            self._truth = self.run_model(self.true_param).ravel()
        return self._truth


    @property
    def cov(self) -> NDArrayFloat:
        if self._cov is None:
            self._cov = np.zeros((self.n_dim, self.n_dim), dtype=self.true_param.dtype)
        return self._cov


    @cov.setter
    def cov(self, c: NDArrayFloat) -> None:
        self._cov = c


    # User can specify anything needed by the model along with the truth and covariance here
    def initialize_model(self) -> None:
        """
        Set all model hyperparameters. Children classes should override
        this method.
        """
        self.A = np.random.normal(loc=0, scale=2, size=(self.n_dim, self.true_param.size))


    # User should calculate the model here. It should be specifically written to run a whole ensemble.
    def run_ensemble(self, param_ens: NDArrayFloat) -> NDArrayFloat:
        """
        Runs a linear model on an ensemble of parameters. Children
        classes should override this method.

        Parameters
        ----------
        param_ens : array
            A 2d array representing an ensemble of model parameters.
        """
        return self.A @ (param_ens)


    # Can run either a single set of parameters or a whole ensemble.
    def run_model(self, parameters: NDArrayFloat) -> NDArrayFloat:
        """
        A more general method of run_ensemble that will work for 
        either an ensemble or single set of parameters.

        Parameters
        ----------
        parameters : array
            A 1d or 2d array of parameters.
        """
        if parameters.ndim == 1:
            parameters = np.expand_dims(parameters, axis=1)
        return self.run_ensemble(parameters)


# Adjust the Model class to include noise
class NoisyLinear(Model):
    """
    A linear model with noise

    Attributes
    ----------
    true_param : array
        A 1d array of target parameters set by the user.
    n_dim : int
        The dimenion of the data space.
    truth : array
        The model response evaluated at the true parameters.
    cov : array
        The covariance of the model response.

    Methods
    -------
    initialize_model()
        Generate or set all hyperparameters for the model.
    run_ensemble(param_ens)
        Run the model on an ensemble of parameters.
    run_model(param)
        Run the model on either a set of parameters or ensemble of
        parameters.
    """
    
    def initialize_model(self) -> None:
        """
        Set all hyperparameters.
        """
        super().initialize_model()
        C = np.random.normal(loc=0, scale=2, size=(self.n_dim, self.n_dim))
        self.cov = 0.05*(C.T.dot(C) + 0.05*np.eye(self.n_dim))


    def run_ensemble(self, param_ens: NDArrayFloat) -> NDArrayFloat:
        """
        Runs a linear model on an ensemble of parameters.

        Parameters
        ----------
        param_ens : array
            A 2d array representing an ensemble of model parameters.
        """
        g = super().run_ensemble(param_ens)
        n_dim, n_ens = g.shape
        eta = np.random.multivariate_normal(np.zeros(n_dim), self.cov, size=n_ens).T
        return g + eta
