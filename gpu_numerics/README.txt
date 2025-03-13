This is a report from my Final Exam project in a GPU Computing course at Caltech.
My work showed that GPUs have promise for Jacobi iteration and multigrid solvers,
but they are tricky to implement. Optimized performance would require going back
and forth between GPU and CPU solvers in V-cycles, which was beyond the scope
of this project.



CPU/GPU Demo: Multigrid Solvers

Emmet M. Cleary


-------------------------------------------------------------------------

1) Outline

2) Compiling

3) Running

4) Testing <-- includes expected results

5) Results

-------------------------------------------------------------------------





-------------------------------------------------------------------------
1) Outline
-------------------------------------------------------------------------

The goal of this project is to solve

-u''(x) = (2 * PI) ^ 2 * cos(2 * PI * x)
u(0) = u(1) = 0

for x in [0,1]. This is accomplished using Jacobi iteration and
multigrid solver.

The code includes functions that are solvers, recursive functions that
call the solvers, and iterative functions that call timers and error
computing functions and print out the results.

SOLVERS:
--------

-- jacobi, jacobiCudaKernel
-- relax, relaxCudaKernel
-- restrict, restrictCudaKernel
-- interpolate, interpolateCudaKernel
-- residual, residualCudaKernel
-- errorCorrectionCudaKernel


RECURSIVE FUNCTIONS:
--------------------

-- vCycle, vCycleCuda
-- fCycle, fCycleCuda


ITERATIVE FUNCTIONS/ERRORS:
---------------------------

-- computeError, computeErrorKernel
-- jacobiSolver, jacobiSolverCuda
-- vCycleSolver, vCycleSolverCuda
-- fCycleSolver, fCycleSolverCuda




-------------------------------------------------------------------------
2) Compiling
-------------------------------------------------------------------------

To compile on Haru,

$ make -f Makefile-haru

I have not tested it on Haru, but the makefiles are based of the
homework, and all of those worked just fine. No external libraries are
needed to compile.


-------------------------------------------------------------------------
3) Running
-------------------------------------------------------------------------

Compiling returns the executable 'multigrid'. It takes 2 inputs:

$ ./multigrid

Usage: multigrid intervals solver [threads]

   intervals: log base 2 of the maximum number of intervals (e.g. 5 yields 2^5)
   solvers: j for Jacobi, v for V-cycle, f for F-cycle
   threads: max number of threads per block

The first input is for the number of intervals in the mesh. An input
of 8 yields a maximum mesh size of 2^8 = 256 intervals (or 257
points), and the code will run the solvers for meshes of 2, 4, 8, ...,
256 intervals.

The solver inputs denote which solver to test.

j: Jacobi
v: V-cycle
f: F-cycle

Threads is optional. Omitting it causes the cpu solver to run, and
including it runs the gpu solvers.


-------------------------------------------------------------------------
4) Testing
-------------------------------------------------------------------------

These solvers filter out different frequencies at different rates. As
the grid size affects the frequencies to be filtered, it is possible
that finer grids will require fewer iterations, hence being faster
than coarser grids (at least when running the CPU demo). IF TIMES DO
NOT SCALE WITH GRID SIZE, THAT IS OKAY.

GPU/CPU Results
---------------

Ideally, the GPU and CPU solvers should give identical results. This
is verified by running different cases. For example,

$ ./multigrid 20 f 512; ./multigrid 20 f

will run the GPU F-cycle for interval numbers ranging from 2 to 2^20
using 512 threads per block, followed by running the same case on
CPU. Iteration counts and errors should be identical.


Convergence
-----------

To verify my solvers are correct, one must check the order of
convergence. Since these solvers are expected to be second order
accurate, we expect that doubling the number of intervals should
decrease the error by a factor of 4.

This is most observable in either the jacobi or v-cycle solver, e.g.

$ ./multigrid 20 v 512

GPU V-cycle solver
L	n	iter	error		time		time/iter
1	2	13	2.93e+00	0.001601	1.23e-04
2	4	14	4.67e-01	0.004180	2.99e-04
3	8	14	1.06e-01	0.006545	4.67e-04
4	16	13	2.59e-02	0.008064	6.20e-04
5	32	13	6.44e-03	0.010551	8.12e-04
6	64	40	1.61e-03	0.039147	9.79e-04
7	128	24	3.99e-04	0.027783	1.16e-03
8	256	20	1.00e-04	0.026380	1.32e-03

Note that it is still correct for the f-cycle, it's just a little less
obvious due to only running a single iteration.

NOTE: Errors do not monotonically decrease, but this is likely due to
the build up of machine errors and me using single precision.



Runtime
-------

As seen above, outputs include both the runtime, but these do not
monotonically increase. This is because these solvers filter out
spectral frequencies which depend on the grid size. Some grids are
just lucky, and filter things out faster (in fewer iterations). The
analysis in the following section will instead focus on the average
runtime (time / iter).



-------------------------------------------------------------------------
5) Results
-------------------------------------------------------------------------


One of the big questions to answer is for what sized problem is GPU
faster than CPU? To study this, times were measured ONLY for the
solvers and error calculations. The solvers did not require moving
data back and forth between CPU and GPU, so even setting up the GPU
solvers were omitted from time measurements.

First, we look at Jacobi iteration solvers. These solvers only use
Jacobi iteration (and are slow at it!), which is a tool utilized by
both cyclic solvers. It is seen that CPU average iteration times a
faster than the GPU, but they also monotonically increase. The GPU
solver is at a roughly constant iteration time average, BECAUSE IT IS
IN PARALLEL! While the GPU times do start increasing around n = 4096,
they outrun the CPU by n = 8192.

GPU Jacobi solver
L	n	iter	error		time		time/iter
1	2	2001	2.93e+00	0.027068	1.35e-05
2	4	2001	4.67e-01	0.027017	1.35e-05
3	8	2001	1.06e-01	0.026997	1.35e-05
4	16	3001	2.59e-02	0.040666	1.36e-05
5	32	5001	6.42e-03	0.067752	1.35e-05
6	64	16001	1.48e-03	0.207825	1.30e-05
7	128	53001	2.32e-04	0.671579	1.27e-05
8	256	174001	2.47e-03	2.230379	1.28e-05
9	512	543001	1.07e-02	7.406110	1.36e-05
10	1024	1589001	4.24e-02	21.706209	1.37e-05
11	2048	3958001	1.74e-01	54.296162	1.37e-05
12	4096	6583001	6.49e-01	102.412842	1.56e-05
13	8192	2304001	1.71e+00	45.119633	1.96e-05  <-- GPU IS FASTER!!!!

CPU Jacobi solver
L	n	iter	error		time		time/iter
1	2	2001	2.93e+00	0.000030	1.50e-08
2	4	2001	4.67e-01	0.000027	1.35e-08
3	8	2001	1.06e-01	0.000043	2.15e-08
4	16	3001	2.59e-02	0.000090	3.00e-08
5	32	5001	6.42e-03	0.000230	4.60e-08
6	64	16001	1.48e-03	0.001390	8.69e-08
7	128	53001	2.32e-04	0.008745	1.65e-07
8	256	174001	2.47e-03	0.045779	2.63e-07
9	512	543001	1.07e-02	0.926822	1.71e-06
10	1024	1589001	4.24e-02	5.441430	3.42e-06
11	2048	3958001	1.74e-01	27.510529	6.95e-06
12	4096	6583001	6.49e-01	93.196320	1.42e-05
13	8192	2304001	1.71e+00	64.968712	2.82e-05 <- CPU BECOMES SLOWER!!!!


We notice the same trend for V-cycles and F-cycles: while the average
iteration times do grow faster on the GPU for these cases than they
did just for Jacobi iteration, ALL THE GPU CASES GROW SLOWER THAN
THEIR CPU COUNTERPARTS.

From this analysis, it is clear that GPU are not as efficient as CPU
for small numbers of intervals. It is also clear that the most
efficient implimentation of multigrid solvers on GPU is to implement
them on both the CPU and GPU together. I did not have the time to
implement that, but it would require careful analysis of time taken to
copy the arrays. The grid where GPU overtakes CPU for speed it not
enough to justify switching to GPU. The CPU must be slower than the
GPU for a sufficiently long time to be worth it. 
