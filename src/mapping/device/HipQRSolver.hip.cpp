#ifdef PRECICE_WITH_HIP

#include "mapping/device/HipQRSolver.hip.hpp"
#include <ginkgo/ginkgo.hpp>
#include <hip/hip_runtime.h>
#include <hip/hip_runtime_api.h>
#include <hipsolver/hipsolver.h>

void computeQRDecompositionHip(const std::shared_ptr<gko::Executor> &exec, GinkgoMatrix *A_Q, GinkgoVector *R)
{
  auto scope_guard = exec->get_scoped_device_id_guard();

  void *dWork{};
  int * devInfo{};

  // Allocating important HIP variables
  hipMalloc((void **) &dWork, sizeof(double));
  hipMalloc((void **) &devInfo, sizeof(int));

  hipsolverDnHandle_t solverHandle;
  hipsolverDnCreate(&solverHandle);
  // NOTE: It's important to transpose since hipsolver assumes column-major memory layout
  // Making a copy since every value will be overridden
  auto A_T = gko::share(gko::matrix::Dense<>::create(exec, gko::dim<2>(A_Q->get_size()[1], A_Q->get_size()[0])));
  A_Q->transpose(A_T);

  // Setting dimensions for solver
  const unsigned int M = A_T->get_size()[1];
  const unsigned int N = A_T->get_size()[0];

  const int lda = max(1, M); // 1 > M ? 1 : M;
  const int k   = max(M, N); // M < N ? M : N;

  int lwork_geqrf = 0;
  int lwork_orgqr = 0;
  int lwork       = 0;

  double *dTau{};
  hipMalloc((void **) &dTau, sizeof(double) * M);

  // Query working space of geqrf and orgqr
  hipsolverStatus_t hipsolverStatus = hipsolverDnDgeqrf_bufferSize(solverHandle, M, N, A_T->get_values(), lda, &lwork_geqrf);
  assert(hipsolverStatus == HIPSOLVER_STATUS_SUCCESS);
  hipsolverStatus = hipsolverDnDorgqr_bufferSize(solverHandle, M, N, k, A_T->get_values(), lda, dTau, &lwork_orgqr);
  assert(hipsolverStatus == HIPSOLVER_STATUS_SUCCESS);
  lwork                   = (lwork_geqrf > lwork_orgqr) ? lwork_geqrf : lwork_orgqr;
  hipError_t hipErrorCode = hipMalloc((void **) &dWork, sizeof(double) * lwork);
  assert(hipSuccess == hipErrorCode);

  void *hWork{};
  // Compute QR factorization
  hipsolverStatus = hipsolverDnDgeqrf(solverHandle, M, N, A_T->get_values(), lda, dTau, (double *) dWork, lwork, devInfo);
  hipErrorCode    = hipDeviceSynchronize();
  assert(hipsolverStatus == HIPSOLVER_STATUS_SUCCESS);
  assert(hipSuccess == hipErrorCode);

  // Copy A_T to R s.t. the upper triangle corresponds to R
  A_T->transpose(R);

  // Compute Q
  hipsolverStatus = hipsolverDnDorgqr(solverHandle, M, N, k, A_T->get_values(), lda, dTau, (double *) dWork, lwork, devInfo);
  hipErrorCode    = hipDeviceSynchronize();
  assert(hipsolverStatus == HIPSOLVER_STATUS_SUCCESS);
  assert(hipSuccess == hipErrorCode);

  A_T->transpose(A_Q);

  hipDeviceSynchronize();

  // Free the utilized memory
  hipFree(dTau);
  hipFree(dWork);
  hipFree(devInfo);
  hipsolverDnDestroy(solverHandle);
}

void solvewithQRDecompositionHip(const std::shared_ptr<gko::Executor> &exec, GinkgoMatrix *U, GinkgoVector *x, GinkgoVector *rhs, GinkgoMatrix *matQ, GinkgoVector *in_vec)
{
  auto scope_guard = exec->get_scoped_device_id_guard();

  hipblasHandle_t handle;
  hipblasStatus_t hipblasStatus = hipblasCreate(&handle);
  assert(hipblasStatus == HIPBLAS_STATUS_SUCCESS);
  double a      = 1;
  double b      = 0;
  hipblasStatus = hipblasDgemv(handle, HIPBLAS_OP_T,
                               matQ->get_size()[0], matQ->get_size()[1],
                               &a,
                               matQ->get_values(), matQ->get_size()[0],
                               in_vec->get_values(), 1,
                               &b,
                               rhs->get_values(), 1);
  assert(hipblasStatus == HIPBLAS_STATUS_SUCCESS);

  hipblasFillMode_t  uplo  = HIPBLAS_FILL_MODE_LOWER;
  hipblasOperation_t trans = HIPBLAS_OP_T;

  // unit triangular = diag = 1
  hipblasDiagType_t diag = HIPBLAS_DIAG_NON_UNIT;
  int               rows = rhs->get_size()[0];
  const int         lda  = max(1, rows);

  hipblasStatus = hipblasDtrsv(handle, uplo,
                               trans, diag,
                               rows, U->get_values(), lda,
                               rhs->get_values(), 1);
  assert(hipblasStatus == HIPBLAS_STATUS_SUCCESS);

  hipDeviceSynchronize();
  *x            = *rhs;
  hipblasStatus = hipblasDestroy(handle);
  assert(hipblasStatus == HIPBLAS_STATUS_SUCCESS);
}
#endif
