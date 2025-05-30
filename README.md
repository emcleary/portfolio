This is a collection of miscellaneous coding projects from my graduate
studies and beyond.

If cloning, I recommend cloning recursively to include all submodules:
```
git clone --recurse-submodules https://github.com/emcleary/portfolio.git
```


## GPU Multigrid Framework (2025)

This is a framework I wrote for a geometric multigrid solver
parallelized with CUDA.

Demonstrated skills:
- C++
- CUDA for GPGPU
- partial differential equations
- iterative solvers
- multigrid methods
- modular code design

## Unstructured Mesh Refinement (2025)

The Unstructured Mesh Refinement (UMR) library generates 2D
unstructured meshes using various Delaunay triangulation and
refinement algorithms.

Demonstrated skills:
- abstract classes (see [sources](https://github.com/emcleary/umr/blob/main/src/sources_shape_interface.hpp), [refinement](https://github.com/emcleary/umr/blob/main/src/refinement_interface.hpp), [parametrics](https://github.com/emcleary/umr/blob/main/src/parametrics_interface.cpp))
- design patterns: [builder](https://github.com/emcleary/umr/blob/main/src/builder.hpp), [factory](https://github.com/emcleary/umr/blob/main/src/quadedge.hpp), [facade](https://github.com/emcleary/umr/blob/main/src/optimizers.hpp)
- data structures: [QuadEdge](https://github.com/emcleary/umr/blob/main/src/quadedge.hpp)
- [command line interface](https://github.com/emcleary/umr/blob/main/src/command_line_interface.hpp)
- external libraries: AlgLib, Boost, Jsoncpp, VTK
- CMake for compiling
- Doxygen for documentation
- Delaunay triangulation

## Game Boy Emulator (2024)

This is not a full release of my emulator, but rather an overview of
my design process.  It includes many excerpts of my code written in
C++, along with discussions of my implementation regarding the
usefulness of specific design patterns and classes.

Demonstrated skills:
- design patterns: factory methods, facade, mediator, state machine
- classes: inheritance, polymorphism
- RAII / resource management
- testing and performance

## PDE Solver with Irregular Geometries (2023)

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

## Optimization (2020)

This is an example test script from my research project implementing
the Ensemble Kalman Inversion algorithm for optimizing model
parameters given noisy data.

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
