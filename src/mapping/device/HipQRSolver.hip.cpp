#ifdef PRECICE_WITH_HIP

#include "mapping/device/HipQRSolver.hip.hpp"
#include <ginkgo/ginkgo.hpp>
#include <hip/hip_runtime.h>
#include <hip/hip_runtime_api.h>
#include <hipsolver.h>

void computeQRDecompositionHip(const int deviceId, const std::shared_ptr<gko::Executor> &exec, gko::matrix::Dense<> *A_Q, gko::matrix::Dense<> *R)
{
  int backupDeviceId{};
  hipGetDevice(&backupDeviceId);
  hipSetDevice(deviceId);

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
  A_Q->transpose(gko::lend(A_T));

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
  A_T->transpose(gko::lend(R));

  // Compute Q
  hipsolverStatus = hipsolverDnDorgqr(solverHandle, M, N, k, A_T->get_values(), lda, dTau, (double *) dWork, lwork, devInfo);
  hipErrorCode    = hipDeviceSynchronize();
  assert(hipsolverStatus == HIPSOLVER_STATUS_SUCCESS);
  assert(hipSuccess == hipErrorCode);

  A_T->transpose(gko::lend(A_Q));

  hipDeviceSynchronize();

  // Free the utilized memory
  hipFree(dTau);
  hipFree(dWork);
  hipFree(devInfo);
  hipsolverDnDestroy(solverHandle);

  // ...and switch back to the GPU used for all coupled solvers
  hipSetDevice(backupDeviceId);
}
#endif
