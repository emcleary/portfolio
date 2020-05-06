#ifndef CUDA_MULTIGRID_CUDA_CUH
#define CUDA_MULTIGRID_CUDA_CUH

void jacobiSolverCuda( const float *f, const float *u, const int n, const unsigned int threadsPerBlock);
void vCycleSolverCuda( const float *f, const float *u, const int n, const unsigned int threadsPerBlock);
void fCycleSolverCuda( const float *f, const float *u, const int n, const unsigned int threadsPerBlock);
float computeErrorCuda( const float *u, const float *v, const int n, float *error, const unsigned int blocks, const unsigned int threadsPerBlock);

#define WEIGHT 0.5

#endif // CUDA_MULTIGRID_CUDA_CUH