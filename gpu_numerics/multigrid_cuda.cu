#include <cstdio>
#include <cuda_runtime.h>
#include "multigrid_cuda.cuh"


//////////////////
// GPU KERNELS //
////////////////

/* Given in problem sets */
__device__
static float atomicMax(float* address, float val)
{
    int* address_as_i = (int*) address;
    int old = *address_as_i, assumed;
    do {
        assumed = old;
        old = ::atomicCAS(address_as_i, assumed,
            __float_as_int(::fmaxf(val, __int_as_float(assumed))));
    } while (assumed != old);
    return __int_as_float(old);
}

/* Adapted from problem sets to find the maximum error using (simple)
   reduction technique. */
__global__
void computeErrorKernel( const float *u, const float *v, const int n, float *error) {

  extern __shared__ float emax[];
  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  // reset error
  if (i == 0) error[0] = 0;

  // reset shared memory
  emax[tid] = 0;

  // Load to shared memory
  while (i < n+1) {
    float val = fabs(u[i] - v[i]);
    emax[tid] = fmax(val, emax[tid]);
    i += gridDim.x * blockDim.x;
  }
  __syncthreads();

  // Find max in block
  /* NOTE: this could be optimized further by unrolling the loop */
  for (unsigned int s = blockDim.x/2; s > 0; s = s / 2) {
    if (tid < s) {
      emax[tid] = fmax(emax[tid], emax[tid+s]);
    }
    __syncthreads();
  }

  // Find max across blocks
  if (tid == 0) atomicMax(error, emax[0]);

}

/* Jacobi solver is the basis of this whole project. Note that several
   different approaches were tested. This approach was taken in all
   subsequent kernels as well, where appropriate, but only the best
   was used. Here, all three were left, if interested. */

__global__
void jacobiCudaKernel( float *v, const float *vprev, const float *f, const int n) {

  extern __shared__ float vprev_sm[];
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  float y;
  float h2 = 1./n/n; // h * h


  // TAKE 1: serial version. Incredibly slow for obvious reasons, but
  // this approach was the first attempt in all my kernels to ensure
  // proper structure of my code.
  /*
  if (i == 0) {
    v[0] = 0;
    for (int j = 1; j<n; j++) {
      y = (vprev[j-1] + vprev[j+1] - f[j]*h2) / 2;
      v[j] = y * (1-WEIGHT) + WEIGHT * vprev[j];
    }
    v[n] = 0;
  }
  */


  // TAKE 2: shared memory. My motivation for this approach was that
  // several indices are needed to compute the desired value. It
  // seemed like a great idea (and was used on Lab 6, problem 1), but
  // the if statements do causes instances of warp
  // divergence. Additionally, reading the endpoints from global
  // memory are inefficient due to memory coalescing.
  /*
  while (i < n+1) {

    // load current data to shared memory
    unsigned int i_sm = threadIdx.x + 1;
    vprev_sm[i_sm] = vprev[i];
    if (threadIdx.x == 0) {
      if (i == 0) {
	vprev_sm[0] = 0;
      } else {
	vprev_sm[0] = vprev[i-1];
      }
    }
    if (threadIdx.x == blockDim.x-1) {
      if (i == n-1) {
	vprev_sm[i_sm+1] = 0;
      } else {
	vprev_sm[i_sm+1] = vprev[i+1];
      }
    }
    __syncthreads();
    
    // Solver
    y = (vprev_sm[i_sm-1] + vprev_sm[i_sm+1] - f[i]*h2) / 2;
    y *= (1-WEIGHT);
    y += WEIGHT * vprev_sm[i_sm];
    v[i] = y;

    // Increment index
    i += blockDim.x * gridDim.x;
  }
  */


  // TAKE 3: Read and write to global memory. Despite global memory
  // coalescing problems, this slightly beat my shared memory version.
  while (i < n) {

    if (i == 0) {
      v[0] = 0;
      v[n] = 0;
    } else {
      y = (vprev[i-1] + vprev[i+1] - f[i]*h2) / 2;
      y *= (1-WEIGHT);
      y += WEIGHT * vprev[i];
      v[i] = y;
      
    }

    // Increment index
    i += blockDim.x * gridDim.x;
  }
  
}

// Computes residual, r = Av - f.
__global__
void residualCudaKernel(float *r, const float *v, const float *f, const int n) {

  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  float n2 = (float) n*n;

  // Load global memory to hide latency. It did not make a big impact,
  // but was worth testing.
  while (i < n) {
    if (i == 0) {
      float v0 = v[0];
      float vn = v[n];
      float f0 = f[0];
      float fn = f[n];
      r[0] = v0 - f0;
      r[n] = vn - fn;
    } else {
      float vim1 = v[i-1];
      float vi   = v[i];
      float vip1 = v[i+1];
      float fi   = f[i];
      r[i]  = (vim1 - 2 * vi + vip1)*n2 - fi;
    }
    i += blockDim.x * gridDim.x;

  }

}

// Coarsens the grid.
__global__
void restrictCudaKernel(float *v2h, float *vh, const int nh) {

  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  
  while (i < nh/2) {
    if (i == 0) {
      v2h[0] = vh[0];
    } else {
      float vh_2im1 = vh[2*i-1];
      float vh_2i   = vh[2*i];
      float vh_2ip1 = vh[2*i+1];
      v2h[i] = (vh_2im1 + 2 * vh_2i + vh_2ip1) / 4;
    }
    i += blockDim.x * gridDim.x;
  }

}

// Refines the grid.
__global__
void interpolateCudaKernel(float *vh, float *v2h, const int n2h) {

  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  while (i < n2h) {
    float v2h_i   = v2h[i];
    float v2h_ip1 = v2h[i+1];
    vh[2*i] = v2h_i;
    vh[2*i+1] = (v2h_i + v2h_ip1) / 2;
    i += blockDim.x * gridDim.x;
  }

}

// Corrects numerical solution with computed errors.
__global__
void errorCorrectionCudaKernel(float *v, const float *e, const int n) {

  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  // I could have unrolled this loop (and many others), but I doubt
  // they were enough of a bottleneck in the code to matter
  // (especially this really simple one!)
  while (i < n+1) {
    v[i] -= e[i];
    i += blockDim.x * gridDim.x;
  }

}

// Iterates over the Jacobi kernel.
void relaxCuda(float *v, float *v0, const float *f, const int nu, const int n,
	       const unsigned int threadsPerBlock) {

  // shared memory parameter
  unsigned int size_sm = (threadsPerBlock + 2) * sizeof(float);
  unsigned int blocks = (int) ceil((float) n / threadsPerBlock);
  
  int it = 0;
  while (it < nu) {

    jacobiCudaKernel<<<blocks, threadsPerBlock, size_sm>>>( v, v0, f, n );

    float *tmp1 = v;
    v  = v0;
    v0 = tmp1;

    it++;
  }

  float *tmp2 = v;
  v  = v0;
  v0 = tmp2;
  
}


/* CYCLES:
 *
 * These are the cyclic solvers that manage changing grid sizes,
 * computing residuals and errors, and using these errors to correct
 * the initial guess (the input solution).
 */
void vCycleCuda(float *v, const float *f, const int L, int l,
		const unsigned int maxThreadsPerBlock) {

  const int n = pow(2, L-l);
  const int nu1 = 2;
  const int nu2 = 1;

  // Compute in this function
  unsigned int threadsPerBlock = min(maxThreadsPerBlock, n);
  unsigned int blocks = (int) ceil((float) n / threadsPerBlock);

  float *vh;
  cudaMalloc( (void **) &vh, (n+1) * sizeof(float) );
  cudaMemset( &vh[0], 0, sizeof(float) );
  cudaMemset( &vh[n], 0, sizeof(float) );

  // Jacobi iteration filtering
  relaxCuda(vh, v, f, nu1, n, threadsPerBlock);
  cudaMemset( &vh[0], 0, sizeof(float) );
  cudaMemset( &vh[n], 0, sizeof(float) );

  if (l < L-1) {

    // Compute residual
    float *rh;
    cudaMalloc( (void **) &rh, (n+1) * sizeof(float) );
    residualCudaKernel<<<blocks, threadsPerBlock>>>(rh, vh, f, n);

    // h -> 2h
    threadsPerBlock = min(n/2, maxThreadsPerBlock);
    blocks = (int) ceil((float) n/2 / threadsPerBlock);
    
    // Restrict residual
    float *r2h;
    cudaMalloc( (void **) &r2h, (n/2+1) * sizeof(float) );
    restrictCudaKernel<<<blocks, threadsPerBlock>>>(r2h, rh, n);
    cudaFree(rh);

    // Next subset of v-cycle: compute error
    float *e2h;
    cudaMalloc( (void **) &e2h, (n/2+1) * sizeof(float) );
    cudaMemset( e2h, 0, (n/2+1) * sizeof(float) );
    vCycleCuda(e2h, r2h, L, l+1, maxThreadsPerBlock);
    cudaFree(r2h);
    
    // Interpolate error
    float *eh;
    cudaMalloc( (void **) &eh, (n+1) * sizeof(float) );
    interpolateCudaKernel<<<blocks, threadsPerBlock>>>(eh, e2h, n/2);
    cudaFree(e2h);

    // 2h -> h
    threadsPerBlock = min(n, maxThreadsPerBlock);
    blocks = (int) ceil((float) n / threadsPerBlock);
    
    // Error correction
    errorCorrectionCudaKernel<<<blocks, threadsPerBlock>>>(vh, eh, n);
    cudaMemset(&vh[0], 0, sizeof(float));
    cudaMemset(&vh[n], 0, sizeof(float));

  }

  // Jacobi iteration filtering
  relaxCuda(v, vh, f, nu2, n, threadsPerBlock);

  cudaFree(vh);
}

void fCycleCuda(float *vh, float *fh, const int L, int l,
		const unsigned int maxThreadsPerBlock) {

  const int n = pow(2, L-l);
  const int nu0 = 1;

  // Compute in this function
  unsigned int blocks;
  unsigned int threadsPerBlock;
  
  if (l < L-1) {

    // h -> 2h
    threadsPerBlock = min(n/2, maxThreadsPerBlock);
    blocks = (int) ceil((float) n/2 / threadsPerBlock);

    // Coarsen input vh
    float *v2h;
    cudaMalloc( (void **) &v2h, (n/2+1) * sizeof(float) );
    restrictCudaKernel<<<blocks, threadsPerBlock>>>(v2h, vh, n);
    cudaMemset(&v2h[0], 0, sizeof(float));
    cudaMemset(&v2h[n/2], 0, sizeof(float));

    // Coarsen forcing function
    float *f2h;
    cudaMalloc( (void **) &f2h, (n/2+1) * sizeof(float) );
    restrictCudaKernel<<<blocks, threadsPerBlock>>>(f2h, fh, n);

    // Coarsen them again, recursively
    fCycleCuda(v2h, f2h, L, l+1, maxThreadsPerBlock);
    cudaFree(f2h);

    // Refine grid for the NEXT V-cycle increment
    interpolateCudaKernel<<<blocks, threadsPerBlock>>>(vh, v2h, n/2);
    cudaFree(v2h);
  }

  // V-cycle!
  int iter = 0;
  while (iter < nu0) {
    vCycleCuda(vh, fh, L, l, maxThreadsPerBlock);
    iter++;
  }
  
}


/* SOLVERS: Called from the host
 *
 * These set up the problem for cuda solvers (e.g. allocate memory),
 * and manage iterations, computing errors, measuring time, and
 * printing results.
 */
void jacobiSolverCuda(const float *f_host, const float *u_host, const int L,
		      const unsigned int threadsPerBlock) {

  const int n = pow(2, L);

  // Numerical solution
  float **v = new float*[2];
  cudaMalloc( (void **) &v[0], (n+1) * sizeof(float) );
  cudaMalloc( (void **) &v[1], (n+1) * sizeof(float) );
  cudaMemset( v[0], 0, (n+1) * sizeof(float) );
  cudaMemset( v[1], 0, (n+1) * sizeof(float) );

  // Forcing function
  float *f;
  cudaMalloc( (void **) &f, (n+1) * sizeof(float) );
  cudaMemcpy( f, f_host, (n+1) * sizeof(float), cudaMemcpyHostToDevice );

  // Exact solution
  float *u;
  cudaMalloc( (void **) &u, (n+1) * sizeof(float) );
  cudaMemcpy( u, u_host, (n+1) * sizeof(float), cudaMemcpyHostToDevice );

  // Error
  float *maxError;
  cudaMalloc( (void **) &maxError, sizeof(float) );
  cudaMemset( maxError, 0, sizeof(float) );

  // Other
  int iter = 0;
  float e_prev = 1;
  float *error = new float[1];
  error[0] = 0;
  unsigned int size_sm = (threadsPerBlock + 2) * sizeof(float);
  unsigned int blocks = (int) ceil( (float) n / threadsPerBlock );

  // Solver
  float tic = clock();
  while (abs(error[0]-e_prev) > 1e-4) {

    float *vprev = v[iter % 2];
    float *vcurr = v[(iter+1) % 2];

    jacobiCudaKernel<<<blocks, threadsPerBlock, size_sm>>>( vcurr, vprev, f, n );

    // Ensure BCs are 0
    cudaMemset( &vcurr[0], 0, sizeof(float) );
    cudaMemset( &vcurr[n], 0, sizeof(float) );

    if (iter % 1000 == 0) {

      // Store previous error
      e_prev = error[0];

      // Compute error
      size_sm = threadsPerBlock * sizeof(float);
      unsigned int blocks = (int) ceil((float) n / threadsPerBlock);
      computeErrorKernel<<<blocks, threadsPerBlock, size_sm>>>( u, vcurr, n, maxError);
      cudaMemcpy( error, maxError, sizeof(float), cudaMemcpyDeviceToHost );
    }

    iter++;
  }
  float toc = clock();
  float time = (toc - tic) / CLOCKS_PER_SEC;
  printf("%d\t%d\t%d\t%6.2e\t%f\t%6.2e\n", L, n, iter, error[0], time, time / iter);

  cudaFree(v[0]);
  cudaFree(v[1]);
  cudaFree(f);
  cudaFree(u);

}

void vCycleSolverCuda(const float *f_host, const float *u_host, const int L,
		      const unsigned int threadsPerBlock) {

  const int n = pow(2, L);

  // Numerical solution
  float *v;
  cudaMalloc( (void **) &v, (n+1) * sizeof(float) );
  cudaMemset( v, 0, (n+1) * sizeof(float) );

  // Forcint function
  float *f;
  cudaMalloc( (void **) &f, (n+1) * sizeof(float) );
  cudaMemcpy( f, f_host, (n+1) * sizeof(float), cudaMemcpyHostToDevice );

  // Exact solution
  float *u;
  cudaMalloc( (void **) &u, (n+1) * sizeof(float) );
  cudaMemcpy( u, u_host, (n+1) * sizeof(float), cudaMemcpyHostToDevice );

  // Error
  float *maxError;
  cudaMalloc( (void **) &maxError, sizeof(float) );
  cudaMemset( maxError, 0, sizeof(float) );

  // Other
  int iter = 0;
  float e_prev = 1;
  float *error = new float[1];
  error[0] = 0;
  unsigned int size_sm = threadsPerBlock * sizeof(float);
  unsigned int blocks = (int) ceil((float) n / threadsPerBlock);

  float tic = clock();
  while (abs(error[0]-e_prev) > 1e-15) {

    vCycleCuda(v, f, L, 0, threadsPerBlock);

    // Ensure BCs are 0
    cudaMemset( &v[0], 0, sizeof(float) );
    cudaMemset( &v[n], 0, sizeof(float) );

    // Store previous error
    e_prev = error[0];

    // Compute error
    computeErrorKernel<<<blocks, threadsPerBlock, size_sm>>>( u, v, n, maxError);
    cudaMemcpy( error, maxError, sizeof(float), cudaMemcpyDeviceToHost );

    iter++;
  }
  float toc = clock();
  float time = (toc - tic) / CLOCKS_PER_SEC;
  printf("%d\t%d\t%d\t%6.2e\t%f\t%6.2e\n", L, n, iter, error[0], time, time / iter);

  cudaFree(v);
  cudaFree(f);
  cudaFree(u);

}


void fCycleSolverCuda(const float *f_host, const float *u_host, const int L,
		      const unsigned int threadsPerBlock) {

  const int n = pow(2, L);

  // Numerical solution
  float *v;
  cudaMalloc( (void **) &v, (n+1) * sizeof(float) );
  cudaMemset( v, 0, (n+1) * sizeof(float) );

  // Forcing function
  float *f;
  cudaMalloc( (void **) &f, (n+1) * sizeof(float) );
  cudaMemcpy( f, f_host, (n+1) * sizeof(float), cudaMemcpyHostToDevice );

  // Exact solution
  float *u;
  cudaMalloc( (void **) &u, (n+1) * sizeof(float) );
  cudaMemcpy( u, u_host, (n+1) * sizeof(float), cudaMemcpyHostToDevice );

  // Error
  float *maxError;
  cudaMalloc( (void **) &maxError, sizeof(float) );
  cudaMemset( maxError, 0, sizeof(float) );

  // Other
  unsigned int size_sm = threadsPerBlock * sizeof(float);
  unsigned int blocks = (int) ceil((float) n / threadsPerBlock);
  
  // Solver
  float tic = clock();
  fCycleCuda(v, f, L, 0, threadsPerBlock);

  // Error
  computeErrorKernel<<<blocks, threadsPerBlock, size_sm>>>( u, v, n, maxError);
  float *error = new float[1];
  cudaMemcpy( error, maxError, sizeof(float), cudaMemcpyDeviceToHost );

  float toc = clock();
  float time = (toc - tic) / CLOCKS_PER_SEC;
  printf("%d\t%d\t%d\t%6.2e\t%f\t%6.2e\n", L, n, 1, error[0], time, time);

  cudaFree(v);
  cudaFree(f);
  cudaFree(u);
  
  
}
