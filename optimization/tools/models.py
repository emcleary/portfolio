from __future__ import annotations
import sys
sys.path.append('../src')
import numpy as np
from argparse import ArgumentParser

from typing import Callable, Optional, Generic
from myTypes import *


class Model(Generic[NDArrayFloat]):
    true_param: NDArrayFloat
    n_variables: int
    _truth: Optional[NDArrayFloat]
    _cov: Optional[NDArrayFloat]

    # n_variables: number of variables
    # parameters: numpy array of true parameters
    def __init__(self, n_variables: int, parameters: NDArrayFloat) -> None:
        self.true_param = parameters
        self.n_variables = n_variables

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
            self._cov = np.zeros((self.n_variables, self.n_variables), dtype=self.true_param.dtype)
        return self._cov


    @cov.setter
    def cov(self, c: NDArrayFloat) -> None:
        self._cov = c


    # User can specify anything needed by the model along with the truth and covariance here
    def initialize_model(self) -> None:
        self.A = np.random.normal(loc=0, scale=2, size=(self.n_variables, self.true_param.size))


    # User should calculate the model here. It should be specifically written to run a whole ensemble.
    def run_ensemble(self, parameters: NDArrayFloat) -> NDArrayFloat:
        return self.A @ (parameters)


    # Can run either a single set of parameters or a whole ensemble.
    def run_model(self, parameters: NDArrayFloat) -> NDArrayFloat:
        if parameters.ndim == 1:
            parameters = np.expand_dims(parameters, axis=1)
        return self.run_ensemble(parameters)


# Adjust the Model class to include noise
class NoisyLinear(Model):
    def initialize_model(self) -> None:
        super().initialize_model()
        C = np.random.normal(loc=0, scale=2, size=(self.n_variables, self.n_variables))
        self.cov = 0.1*(C.T.dot(C) + 0.05*np.eye(self.n_variables))


    def run_ensemble(self, parameters: NDArrayFloat) -> NDArrayFloat:
        g = super().run_ensemble(parameters)
        n_variables, n_particles = g.shape
        eta = np.random.multivariate_normal(np.zeros(n_variables), self.cov, size=n_particles).T
        return g + eta
