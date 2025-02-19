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

void computeQRDecompositionCuda(const std::shared_ptr<gko::Executor> &exec, GinkgoMatrix *A_Q, GinkgoVector *R)
{
  auto scope_guard = exec->get_scoped_device_id_guard();

  void *dWork{};
  int * devInfo{};

  // Allocating important CUDA variables
  cudaError_t cudaErrorCode = cudaMalloc((void **) &dWork, sizeof(double));
  assert(cudaSuccess == cudaErrorCode);
  cudaErrorCode = cudaMalloc((void **) &devInfo, sizeof(int));
  assert(cudaSuccess == cudaErrorCode);

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
  cudaErrorCode = cudaMalloc((void **) &dTau, sizeof(double) * M);
  assert(cudaSuccess == cudaErrorCode);

  precice::profiling::Event calculateQRDecompEvent{"calculateQRDecomp"};

  // Query working space of geqrf and orgqr
  cusolverStatus_t cusolverStatus = cusolverDnXgeqrf_bufferSize(solverHandle, nullptr, M, N, CUDA_R_64F, A_T->get_values(), lda, CUDA_R_64F, dTau, CUDA_R_64F, &dLwork_geqrf, &hLwork_geqrf);
  // PRECICE_ASSERTs collide with cuda for some (non-extensively investigated) reason
  assert(cusolverStatus == CUSOLVER_STATUS_SUCCESS);
  cusolverStatus = cusolverDnDorgqr_bufferSize(solverHandle, M, N, k, A_T->get_values(), lda, dTau, (int *) &dLwork_orgqr);
  assert(cusolverStatus == CUSOLVER_STATUS_SUCCESS);
  dLwork                    = (dLwork_geqrf > dLwork_orgqr) ? dLwork_geqrf : dLwork_orgqr;
  cudaErrorCode = cudaMalloc((void **) &dWork, sizeof(double) * dLwork);
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

  cudaErrorCode = cudaDeviceSynchronize();
  assert(cudaSuccess == cudaErrorCode);
  calculateQRDecompEvent.stop();

  // Free the utilized memory
  cudaErrorCode = cudaFree(dTau);
  assert(cudaSuccess == cudaErrorCode);
  cudaErrorCode = cudaFree(dWork);
  assert(cudaSuccess == cudaErrorCode);
  cudaErrorCode = cudaFree(devInfo);
  assert(cudaSuccess == cudaErrorCode);
  cusolverStatus = cusolverDnDestroy(solverHandle);
  assert(cusolverStatus == CUSOLVER_STATUS_SUCCESS);
}

void solvewithQRDecompositionCuda(const std::shared_ptr<gko::Executor> &exec, GinkgoMatrix *U, GinkgoVector *x, GinkgoVector *rhs, GinkgoMatrix *matQ, GinkgoVector *in_vec)
{
  auto scope_guard = exec->get_scoped_device_id_guard();

  cublasHandle_t handle;
  cublasStatus_t cublasStatus = cublasCreate(&handle);
  assert(cublasStatus == CUBLAS_STATUS_SUCCESS);
  double a     = 1;
  double b     = 0;
  cublasStatus = cublasDgemv(handle, CUBLAS_OP_T,
                             matQ->get_size()[0], matQ->get_size()[1],
                             &a,
                             matQ->get_values(), matQ->get_size()[0],
                             in_vec->get_values(), 1,
                             &b,
                             rhs->get_values(), 1);
  assert(cublasStatus == CUBLAS_STATUS_SUCCESS);

  cublasFillMode_t  uplo  = CUBLAS_FILL_MODE_LOWER;
  cublasOperation_t trans = CUBLAS_OP_T;

  // unit triangular = diag = 1
  cublasDiagType_t diag    = CUBLAS_DIAG_NON_UNIT;
  int              rows    = rhs->get_size()[0];
  const int        lda     = max(1, rows);

  cublasStatus = cublasDtrsv(handle, uplo,
                             trans, diag,
                             rows, U->get_values(), lda,
                             rhs->get_values(), 1);
  assert(cublasStatus == CUBLAS_STATUS_SUCCESS);

  // In case we refactor the code in the future to make use of
  // dtrsm instead of dtrsv (processing vector data as a whole),
  // the following holds
  // double           alpha   = 1.0;
  // int              columns = 1;
  // const int        ldb     = max(1, rows);
  // cublasSideMode_t  side  = CUBLAS_SIDE_LEFT;
  // cublasStatus = cublasDtrsm(handle,
  //                            side,
  //                            uplo,
  //                            trans,
  //                            diag,
  //                            rows,
  //                            columns,
  //                            &alpha,
  //                            U->get_values(), lda,
  //                            rhs->get_values(),
  //                            ldb);

  cudaDeviceSynchronize();
  *x           = *rhs;
  cublasStatus = cublasDestroy(handle);
  assert(cublasStatus == CUBLAS_STATUS_SUCCESS);
}
#endif
