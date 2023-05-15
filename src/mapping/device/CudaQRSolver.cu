#ifdef PRECICE_WITH_CUDA

#include <cublas_v2.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <cusolverDn.h>
#include <ginkgo/ginkgo.hpp>
#include "device_launch_parameters.h"
#include "mapping/device/CudaQRSolver.cuh"
#include "profiling/Event.hpp"
#include "profiling/EventUtils.hpp"

void computeQRDecompositionCuda(const int deviceId, const std::shared_ptr<gko::Executor> &exec, gko::matrix::Dense<> *A_Q, gko::matrix::Dense<> *R)
{
  int backupDeviceId{};
  cudaGetDevice(&backupDeviceId);
  cudaSetDevice(deviceId);

  void *dWork{};
  int * devInfo{};

  // Allocating important CUDA variables
  cudaMalloc((void **) &dWork, sizeof(double));
  cudaMalloc((void **) &devInfo, sizeof(int));

  cusolverDnHandle_t solverHandle;
  cusolverDnCreate(&solverHandle);
  // NOTE: It's important to transpose since cuSolver assumes column-major memory layout
  // Making a copy since every value will be overridden
  auto A_T = gko::share(gko::matrix::Dense<>::create(exec, gko::dim<2>(A_Q->get_size()[1], A_Q->get_size()[0])));
  A_Q->transpose(gko::lend(A_T));

  // Setting dimensions for solver
  const unsigned int M = A_T->get_size()[1];
  const unsigned int N = A_T->get_size()[0];

  const int lda = max(1, M);
  const int k   = min(M, N);

  size_t dLwork_geqrf = 0;
  size_t dLwork_orgqr = 0;
  size_t dLwork       = 0;

  size_t hLwork_geqrf = 0;
  size_t hLwork       = 0;

  double *dTau{};
  cudaMalloc((void **) &dTau, sizeof(double) * M);

  precice::profiling::Event calculateQRDecompEvent{"calculateQRDecomp"};

  // Query working space of geqrf and orgqr
  cusolverStatus_t cusolverStatus = cusolverDnXgeqrf_bufferSize(solverHandle, nullptr, M, N, CUDA_R_64F, A_T->get_values(), lda, CUDA_R_64F, dTau, CUDA_R_64F, &dLwork_geqrf, &hLwork_geqrf);
  // PRECICE_ASSERTs collide with cuda for some (non-extensively investigated) reason
  assert(cusolverStatus == CUSOLVER_STATUS_SUCCESS);
  cusolverStatus = cusolverDnDorgqr_bufferSize(solverHandle, M, N, k, A_T->get_values(), lda, dTau, (int *) &dLwork_orgqr);
  assert(cusolverStatus == CUSOLVER_STATUS_SUCCESS);
  dLwork                    = (dLwork_geqrf > dLwork_orgqr) ? dLwork_geqrf : dLwork_orgqr;
  cudaError_t cudaErrorCode = cudaMalloc((void **) &dWork, sizeof(double) * dLwork);
  assert(cudaSuccess == cudaErrorCode);

  void *hWork{};
  // Compute QR factorization
  cusolverStatus = cusolverDnXgeqrf(solverHandle, nullptr, M, N, CUDA_R_64F, A_T->get_values(), lda, CUDA_R_64F, dTau, CUDA_R_64F, dWork, dLwork, hWork, hLwork, devInfo);
  cudaErrorCode  = cudaDeviceSynchronize();
  assert(cusolverStatus == CUSOLVER_STATUS_SUCCESS);
  assert(cudaSuccess == cudaErrorCode);

  // Copy A_T to R s.t. the upper triangle corresponds to R
  A_T->transpose(gko::lend(R));

  // Compute Q
  cusolverStatus = cusolverDnDorgqr(solverHandle, M, N, k, A_T->get_values(), lda, dTau, (double *) dWork, dLwork, devInfo);
  cudaErrorCode  = cudaDeviceSynchronize();
  assert(cusolverStatus == CUSOLVER_STATUS_SUCCESS);
  assert(cudaSuccess == cudaErrorCode);

  A_T->transpose(gko::lend(A_Q));

  cudaDeviceSynchronize();
  calculateQRDecompEvent.stop();

  // Free the utilized memory
  cudaFree(dTau);
  cudaFree(dWork);
  cudaFree(devInfo);
  cusolverDnDestroy(solverHandle);

  // ...and switch back to the GPU used for all coupled solvers
  cudaSetDevice(backupDeviceId);
}
#endif
