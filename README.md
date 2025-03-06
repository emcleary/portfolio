This is a collection of miscellaneous coding projects from my graduate
studies and beyond.

## Unstructured Mesh Refinement (2025)

This library generates 2D unstructured meshes using various Delaunay triangulation and refinement algorithms.

Demonstrated skill:
- abstract classes
- design patterns: builder, factory, facade
- command line interface
- external libraries: AlgLib, Boost, Jsoncpp, VTK
- CMake for compiling
- Doxygen for documentation
- Delaunay triangulation

## PDE Solver with Irregular Geometries (2022, 2023)

This implements a finite difference scheme to solve a PDE on an
irregular geometry. The code can work with different geometries simply
by defining a new child class defining the geometry.

Demonstrated skills (C++ version):
- class inheritance
- builder design pattern
- OpenMP for parallelizing bottleneck calculations
- external libraries Boost and GSL
- CMake and Make for compiling
- Doxygen for documentation
- partial differential equations
- finite difference schemes

Demonstrated skills (Python version):
- class inheritance
- NumPy vectorization for efficient calculations
- external libraries matplotlib and scipy
- breath first search algorithm
- partial differential equations
- finite difference schemes

## Optimization (2020, 2022)

This is an example test script from my research project implementing
the Ensemble Kalman Inversion algorithm for optimizing model
parameters given noisy data. It was updated in 2022 to include
docstrings and type hinting for documentation purposes.

Demonstrated skills:
- class inheritance
- NumPy vectorization for efficient calculations
- pytest for unit testing
- docstrings and type-hinting for documentation
- optimization algorithms (derivative-free/ensemble based)

## Asynchrony Tolerant Schemes (2017)

Asynchrony tolerant schemes are numerical finite difference schemes
designed to reduce latency between sending/receiving ghost cells from
neighboring cores. This excerpt is a set of functions I wrote in
Fortran to use these schemes to solve 1D transport equations for
reacting flows. The functions call other in-house codes/libraries
(excluded) for the numerical schemes, calculations of transport
properties, and calculations of reactions rates.

This work served as part of a preliminary study for [this journal
article](https://www.sciencedirect.com/science/article/abs/pii/S0021999123000013).

Demonstrated skills:
- Fortran (f77)
- finite difference schemes
- partial differential equations
- fluid mechanics and combustion problems

## GPU Numerics (2017)

This was a final project I wrote for a GPU programming course. The
project solves an ODE using Jacobi iteration and multigrid solvers on
both CPU and GPU using CUDA.

Demonstrated skills:
- C
- CUDA for GPGPU
- partial differential equations
- iterative solvers
- multigrid methods
