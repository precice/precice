#pragma once
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <cusolverDn.h>
#include <ginkgo/ginkgo.hpp>

using GinkgoMatrix = gko::matrix::Dense<>;

class QRSolver {
  public:
   QRSolver(const int deviceId = 0);
   void computeQR(const std::shared_ptr<gko::Executor> &exec, GinkgoMatrix *A_Q, GinkgoMatrix *R);
  ~QRSolver();
private:
  // Handles for low-level CUDA libraries
  cusolverDnHandle_t solverHandle;
  cusolverStatus_t   cusolverStatus = CUSOLVER_STATUS_SUCCESS;
  cudaError_t        cudaErrorCode  = cudaSuccess;

  // Important variables which track the state of the solver routines
  double *dTau    = nullptr;
  void   *dWork   = nullptr;
  void   *hWork   = nullptr;
  int    *devInfo = nullptr;

  int cudaBackupDeviceId = 0;
};