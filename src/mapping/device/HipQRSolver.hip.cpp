#ifdef PRECICE_WITH_HIP

#include "mapping/device/HipQRSolver.hip.hpp"
#include <ginkgo/ginkgo.hpp>
#include <hip/hip_runtime.h>
#include <hip/hip_runtime_api.h>
#include <hipsolver.h>


HipQRSolver::HipQRSolver(const int deviceId)
{
  hipGetDevice(&backupDeviceId);
  hipSetDevice(deviceId);

  // Allocating important HIP variables
  hipMalloc((void **) &dWork, sizeof(double));
  hipMalloc((void **) &devInfo, sizeof(int));
  hipMalloc((void **) &dTau, sizeof(double));

  hipsolverDnCreate(&solverHandle);
}

HipQRSolver::~HipQRSolver()
{
  hipFree(dTau);
  hipFree(dWork);
  hipFree(devInfo);

  hipsolverDnDestroy(solverHandle);
}

void HipQRSolver::computeQR(const std::shared_ptr<gko::Executor> &exec, GinkgoMatrix *A_Q, GinkgoMatrix *R)
{
  // NOTE: It's important to transpose since hipsolver assumes column-major memory layout
  // Making a copy since every value will be overridden
  auto A_T = gko::share(GinkgoMatrix::create(exec, gko::dim<2>(A_Q->get_size()[1], A_Q->get_size()[0])));
  A_Q->transpose(gko::lend(A_T));

  // Setting dimensions for solver
  const unsigned int M = A_T->get_size()[1];
  const unsigned int N = A_T->get_size()[0];

  const int lda = 1 > M ? 1 : M;
  const int k   = M < N ? M : N;

  int lwork_geqrf = 0;
  int lwork_orgqr = 0;
  int lwork       = 0;

  // Query working space of geqrf and orgqr
  hipsolverStatus = hipsolverDnDgeqrf_bufferSize(solverHandle, M, N, A_T->get_values(), lda, &lwork_geqrf);
  assert(hipsolverStatus == HIPSOLVER_STATUS_SUCCESS);
  hipsolverStatus = hipsolverDnDorgqr_bufferSize(solverHandle, M, N, k, A_T->get_values(), lda, dTau, &lwork_orgqr);
  assert(hipsolverStatus == HIPSOLVER_STATUS_SUCCESS);
  lwork         = (lwork_geqrf > lwork_orgqr) ? lwork_geqrf : lwork_orgqr;
  hipErrorCode = hipMalloc((void **) &dWork, sizeof(double) * lwork);
  assert(hipSuccess == hipErrorCode);

  // Compute QR factorization
  hipsolverStatus = hipsolverDnDgeqrf(solverHandle, M, N, A_T->get_values(), lda, dTau, (double *) dWork, lwork, devInfo);
  hipErrorCode    = hipDeviceSynchronize();
  assert(hipsolverStatus == HIPSOLVER_STATUS_SUCCESS);
  assert(hipSuccess == hipErrorCode);

  // Copy A_T to R s.t. the upper triangle corresponds to R
  A_T->transpose(gko::lend(R));

  // Compute Q
  hipsolverStatus = hipsolverDnDorgqr(solverHandle, M, N, k, A_T->get_values(), lda, dTau, (double *) dWork, lwork, devInfo);
  hipErrorCode   = hipDeviceSynchronize();
  assert(hipsolverStatus == HIPSOLVER_STATUS_SUCCESS);
  assert(hipSuccess == hipErrorCode);

  A_T->transpose(gko::lend(A_Q));

  hipDeviceSynchronize();
  hipSetDevice(backupDeviceId); // Switch back to the GPU used for all coupled solvers
}
#endif
