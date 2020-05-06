This is a collection of miscellaneous coding projects from my graduate
studies.


**Optimization**

This is an example test script from my research project, implementing
the Ensemble Kalman Inversion algorithm for optimizing noisy data.

**Asynchrony Tolerant Schemes**

Asynchrony tolerant schemes are numerical finite difference schemes
designed to reduce latency between sending/receiving ghost cells from
neighboring cores. This excerpt is a set of functions I wrote in
Fortran to use these schemes to solve 1D transport equations for
reacting flows. The functions call other in-house codes/libraries
(excluded) for the numerical schemes, calculations of transport
properties, and calculations of reactions rates.

**GPU Numerics**

This was a final project I wrote for a GPU programming course. The
project solves the ODE

```
-u''(x) = (2 * PI) ^ 2 * cos(2 * PI * x)

u(0) = u(1) = 0
```

using Jacobi iteration and multigrid solvers on both CPU and GPU using
CUDA.

**Wrappers**

This was a group project for a course on scientific
computing. Chemical reactions have what is knowns as a progress
variable to track how close the reaction is to completion. This code
check different progress variables by finding the best
combination. The code was written in C++ with SWIG wrappers to Python
for simpler unit testing.

