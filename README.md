This is a collection of miscellaneous coding projects from my graduate
studies and since then. 


## PDE Solver with Irregular Geometries (2022)

This implements a finite difference scheme to solve a PDE on an
irregular geometry. The code is robust and is easy to test with
different geometries.

## Optimization (2020, 2022)

This is an example test script from my research project implementing
the Ensemble Kalman Inversion algorithm for optimizing model
parameters given noisy data. It was updated in 2022 to include
docstrings and type hinting for documentation purposes.

## Asynchrony Tolerant Schemes (2017)

Asynchrony tolerant schemes are numerical finite difference schemes
designed to reduce latency between sending/receiving ghost cells from
neighboring cores. This excerpt is a set of functions I wrote in
Fortran to use these schemes to solve 1D transport equations for
reacting flows. The functions call other in-house codes/libraries
(excluded) for the numerical schemes, calculations of transport
properties, and calculations of reactions rates.

## GPU Numerics (2017)

This was a final project I wrote for a GPU programming course. The
project solves an ODE using Jacobi iteration and multigrid solvers on
both CPU and GPU using CUDA.

## Wrappers (2015)

This was a group project for a course on scientific computing. The
code was primarily written in C++ for speed and called using Python
code for unit testing by writing wrappers with SWIG. The project
itself was used to study chemical mechanism by tracking their progress.
