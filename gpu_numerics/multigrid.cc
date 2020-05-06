#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <fstream>
#include <string.h>

#include <cuda_runtime.h>
#include <algorithm>

#include "multigrid_cuda.cuh"
#include "ta_utilities.hpp"


using namespace std;


// This is a single jacobi iteration. The Cuda version will have
// threads load v to shared memory, use the shared memory to compute
// the updates, then rescale the global memory with the weight w and
// add the updated value to global memory.
void jacobi(float *vnew, const float *v, const float *f,
	    const int n) {

  float h2 = 1./n/n; // h * h
  for (int i = 1; i < n; i++) {
    vnew[i] = (v[i-1] + v[i+1] - f[i]*h2) / 2; // "updated" value
    vnew[i] = vnew[i] * (1-WEIGHT) + WEIGHT * v[i]; // store in global memory
  }

}

// This function calls jacobi for nu iterations. The GPU version will
// just be a simple wrapper to call the jacobi Cuda kernel. Since
// jacobi will use shared memory, this relax/jacobi wrapper will not
// need to create new arrays and copy values.
void relax(float *vnew, float *v, const float *f, const int nu, const int n) {

  int it = 0;
  while (it < nu) {
    jacobi(vnew, v, f, n);
    swap(vnew, v);
    it++;
  }
  swap(vnew, v);

}


// This function coarsens the mesh, reducing the mesh from n intervals
// to n/2. Again, it will require loading from global memory to shared
// memory, computing a value, and storing back in global memory but
// ONLY to half the elements of the array. 
void restrict(float *v2h, const float *vh, const int nh) {

  int n2h = nh/2;
  for (int i = 1; i < n2h; i++){
    v2h[i] = (vh[2*i-1] + 2 * vh[2*i] + vh[2*i+1]) / 4;
  }

  v2h[0] = vh[0];
  v2h[n2h] = vh[nh];
    
}

// This function refines the mesh, doubling the number of
// intervals. Again, I'll use shared memory, but it will write to
// twice the amount of global memory as it reads.
void interpolate(float *vh, const float *v2h, const int n2h) {
  
  for (int i = 0; i < n2h; i++) {
    vh[2*i] = v2h[i];
    vh[2*i+1] = (v2h[i] + v2h[i+1]) / 2;
  }
  
}

// This function computed the residual, i.e. the difference between Av
// and f (left hand side and right hand side of the system of
// equations). Again, this will require shared memory.
void residual(float *r, const float *v, const float *f, const int n) {

  float h = 1./n;
  float Ajm1 =  1./h/h;
  float Aj   = -2./h/h;
  float Ajp1 =  1./h/h;
  r[0] = v[0] - f[0];
  for (int i = 1; i < n; i++) {
    r[i] = v[i-1] * Ajm1 + v[i] * Aj + v[i+1] * Ajp1 - f[i];
  }
  r[n] = v[n] - f[n];

}

// This recursive function calls the above functions: relax, residual,
// restrict, and interpolate.
//
// 1) "relax" by executing nu1 jacobi iterations to compute v
// 2) if coarsest mesh, skip to 6
// 3) compute the residual with this new v
// 4) coarsen the mesh with restrict
// 5) call vCycle on coarsen mesh, return to 1
// 6) "relax" by executing nu2 jacobi iterations to compute errors (v)
// 7) refine the mesh with interpolate to compute error (v)
// 8) correct v with computed error, return to 6 if not finest mesh
// 7) "relax" by executing nu2 jacobi iterations to compute v (no longer error!)
//
// For GPU, this will simply be a wrapper that calls the relevant GPU
// kernels.
void vCycle(float *v, const float *f, const int L, int l) {

  int n = pow(2, L-l);
  int nu1 = 2;
  int nu2 = 1;
  
  float *vh  = new float[n+1];
  vh[0] = 0; vh[n] = 0;
  relax(vh, v, f, nu1, n);
  if (l < L-1) {

    // residual
    float *rh = new float[n+1];
    residual(rh, vh, f, n);

    // Restrict residual
    float *r2h = new float[n/2+1];
    restrict(r2h, rh, n);
    delete rh;

    // next subset of v-cycle: compute error
    float *e2h = new float[n/2+1]();
    vCycle(e2h, r2h, L, l+1);
    delete r2h;
    
    // Interpolate error
    float *eh  = new float[n+1];
    interpolate(eh, e2h, n/2);
    delete e2h;

    // Error correction
    for (int i = 0; i < n; i++) {
      vh[i] = vh[i] - eh[i];
    }
    delete eh;
  }

  relax(v, vh, f, nu2, n);
  delete vh;

}

// This function coarsens the mesh to the finest mesh, then computes v
// using v-cycles. When a v-cycle is complete, it refines the mesh by
// 1 step, the computes another v-cycle, again and again, until the
// finest mesh is reached. Then, it computes one last v-cycle, and
// then is finished.
//
// This solver might be too fast for my GPU codes to be worth
// it. We'll see! Provided I get my GPU v-cycles working, it should be
// quick and easy enough to write this up and see.
void fCycle(float *vh, const float *fh, const int L, int l) {

  int n = pow(2, L-l);
  int nu0 = 1;

  if (l < L-1) {
    float *v2h = new float[n/2+1];
    restrict(v2h, vh, n);

    float *f2h = new float[n/2+1];
    restrict(f2h, fh, n);

    fCycle(v2h, f2h, L, l+1);
    delete f2h;

    interpolate(vh, v2h, n/2);
    delete v2h;
  }

  int iter = 0;
  while (iter < nu0) {
    vCycle(vh, fh, L, l);
    iter++;
  }
  
}

// Computes the maximum error. On GPU, this will be parallelized using
// a simple reduction method.
float computeError(const float *v, const  float *u, const int n) {
  float e = 0;
  for (int i = 0; i < n+1; i++) {
    e = max(abs(v[i]-u[i]), e);
  }
  return e;
}

// This function calls the jacobi iterations until the difference
// between steps of the iterations is below 1e-6. This function
// computes the runtime and the error. For GPU, the code will be
// nearly identical to this version, except it will call the GPU
// wrappers for the jacobi and computeError kernels.
void jacobiSolver(const float *f, const float *u, const int L) {

  const int n = pow(2, L);
  float *v = new float[n+1]();
  float *v0 = new float[n+1]();
  
  int iter = 0;
  float e_prev = 1;
  float e_curr = 0;
  float tic = clock();
  while (abs(e_prev-e_curr) > 1e-4) {
    jacobi(v, v0, f, n);
    if (iter % 1000 == 0) {
      e_prev = e_curr;
      e_curr = computeError(u, v, n);
    }
    swap(v, v0);
    iter++;
  }
  float toc = clock();
  float time = (toc-tic) / CLOCKS_PER_SEC;
  printf("%d\t%d\t%d\t%6.2e\t%f\t%6.2e\n", L, n, iter, e_curr, time, time / iter);

  delete v;
  delete v0;

}


// This function iteratively calls vCycle until the error is 1e-15, as
// well as measures the runtime. It will also be almost identical to
// its GPU counterpart, just calling GPU wrappers for vCycle and
// computeError.
void vCycleSolver(const float *f, const float *u, const int L) {

  const int n = pow(2, L);
  float *v = new float[n+1](); // Initial guess set to 0; could have been something else, like rand;

  int iter = 0;
  float e_curr = 0;
  float e_prev = 1;
  float tic = clock();
  while (abs(e_prev-e_curr) > 1e-15) {
    vCycle(v, f, L, 0);
    e_prev = e_curr;
    e_curr = computeError(u, v, n);
    iter++;
  }
  float toc = clock();
  float time = (toc-tic) / CLOCKS_PER_SEC;
  printf("%d\t%d\t%d\t%6.2e\t%f\t%6.2e\n", L, n, iter, e_curr, time, time / iter);
  
  delete v;

}

// This function only calls fCycle once, because that's all that's
// needed! It also measures the runtime, and it differs from it's GPU
// counterpart only by fCycle and computeError.
void fCycleSolver(const float *f, const float *u, const int L) {

  const int n = pow(2, L);
  float *v = new float[n+1](); // Initial guess set to 0; could have been something else, like rand;
  float tic = clock();
  fCycle(v, f, L, 0);
  float error = computeError(u, v, n);

  float toc = clock();
  float time = (toc-tic) / CLOCKS_PER_SEC;
  printf("%d\t%d\t%d\t%6.2e\t%f\t%6.2e\n", L, n, 1, error, time, time);
  
  delete v;

}


int main(int argc, char* argv[]) {

  // These functions allow you to select the least utilized GPU
  // on your system as well as enforce a time limit on program execution.
  // Please leave these enabled as a courtesy to your fellow classmates
  // if you are using a shared computer. You may ignore or remove these
  // functions if you are running on your local machine.

  TA_Utilities::select_least_utilized_GPU();
  int max_time_allowed_in_seconds = 300;
  TA_Utilities::enforce_time_limit(max_time_allowed_in_seconds);

  // Usage
  if (argc == 3) {
    printf("\nTesting CPU solver with max n = %d\n", int(pow(2, atoi(argv[1]))));
  } else if (argc == 4) {
    printf("\nTesting GPU solver with max n = %d\n", int(pow(2, atoi(argv[1]))));
  } else {
    printf("\nUsage: multigrid intervals solver [threads]\n\n");    
    printf("   intevals: log base 2 of the maximum number of intervals (e.g. 8 yields 2^8 = 256)\n");
    printf("   solvers: j for Jacobi, v for V-cycle, f for F-cycle\n");
    printf("   threads: max number of threads per block\n\n");
    exit(-1);
  }

  int maxL = atoi(argv[1]);
  int *L = new int[maxL];
  for (int i = 0; i < maxL; i++) {
    L[i] = i+1;
  }

  string cases = argv[2];
  unsigned int threadsPerBlock = 0;
  if (argc == 4) {
    threadsPerBlock = atoi(argv[3]);
  }

  // NOTE: The following will look innefficient to you because it has
  // to recompute the forcing function and the exact solution at each
  // grid size for each case. This is true, but doing it this way made
  // my printouts far cleaner. Furthermore, I am only measuring the
  // runtime of the solvers, because that is all that matters (and
  // copying to the GPU, of course!).
  
  // Loops through grid sizes
  for (int nn = 0; nn < maxL; nn++) {
    
    int n = pow(2, L[nn]);
    float h = 1./n;
    
    // Compute RHS
    float *f = new float[n+1];
    for (int i = 0; i < n+1; i++) {
      f[i] = -4 * M_PI * M_PI * cos ( 2 * M_PI * i * h );
    }
    
    // Compute exact solution
    float *exact = new float[n+1];
    for (int i = 0; i < n+1; i++) {
      exact[i] = cos ( 2 * M_PI * i * h ) - 1;
    }

    // CPU cases
    if (threadsPerBlock == 0) {
      if (cases == "j") {
	if (L[nn] == 1) {
	  printf("\nCPU Jacobi solver\n");
	  printf("L\tn\titer\terror\t\ttime\t\ttime/iter\n");
	}
	
	// Jacobi iteration
	jacobiSolver(f, exact, L[nn]);
	
      } else if (cases == "v") {
	if (L[nn] == 1) {
	  printf("\nCPU V-cycle solver\n");
	  printf("L\tn\titer\terror\t\ttime\t\ttime/iter\n");
	}	  
	
	// CPU V-Cycles
	vCycleSolver(f, exact, L[nn]);
	
      } else if (cases == "f") {
	if (L[nn] == 1) {
	  printf("\nCPU F-cycle solver\n");
	  printf("L\tn\titer\terror\t\ttime\t\ttime/iter\n");
	}	  
	
	// CPU F-Cycles
	fCycleSolver(f, exact, L[nn]);
	
      } else {
	printf("Wrong case doesn't exist!");
	exit(1);
      }

      // GPU cases
    } else {
      
      if (cases == "j") {
	if (L[nn] == 1) {
	  printf("\nGPU Jacobi solver\n");
	  printf("L\tn\titer\terror\t\ttime\t\ttime/iter\n");
	}
	
	// GPU Jacobi iteration
	jacobiSolverCuda(f, exact, L[nn], threadsPerBlock);
	
      } else if (cases == "v") {
	if (L[nn] == 1) {
	  printf("\nGPU V-cycle solver\n");
	  printf("L\tn\titer\terror\t\ttime\t\ttime/iter\n");
	}	  
	
	// GPU V-Cycles
	vCycleSolverCuda(f, exact, L[nn], threadsPerBlock);
	
      } else if (cases == "f") {
	if (L[nn] == 1) {
	  printf("\nGPU F-cycle solver\n");
	  printf("L\tn\titer\terror\t\ttime\t\ttime/iter\n");
	}	  
	
	// GPU F-Cycles
	fCycleSolverCuda(f, exact, L[nn], threadsPerBlock);
	
      } else {
	printf("Wrong case doesn't exist!");
	exit(1);
      }
      
    }
    delete f;
    delete exact;
  
  }
  printf("\n");
  return 0;
}
  


