# Ensemble Kalman Inversion

Tools for Ensemble Kalman Inversion (EKI), a derivative free
optimization algorithm useful for working with noisy data.

This algorithm is particularly useful for working with computationally
expensive models, as model evaluations for each iteration are
embarassingly parallelizable.

# Features

This tool has a few handy features:
* It has a tool for experimenting with linear models from the commandline (see `tools/run_cmd_line.py`).
* It has a template Model class that can be used for building additional, simple models for testing.
* EKI can be run either with instances of a Model class or by streaming data (see `examples`).
* The code includes extensive documentation and type hints.

# References

Iglesias et al. Ensemble Kalman methods for inverse problems. Inverse
Problems 29. 045001 (2013)
