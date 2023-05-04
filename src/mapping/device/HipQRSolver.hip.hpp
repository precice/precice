#ifdef PRECICE_WITH_HIP
#pragma once

#include <ginkgo/ginkgo.hpp>
#include "mapping/QRSolver.hpp"
#ifdef __HIPCC__
#ifndef __HIP_PLATFORM_AMD__
#define __HIP_PLATFORM_AMD__
#endif
#include <hip/hip_runtime.h>
#include <hip/hip_runtime_api.h>
#include <hipsolver.h>
#endif
class HipQRSolver : public QRSolver {
public:
  HipQRSolver(const int deviceId = 0);
  void computeQR(const std::shared_ptr<gko::Executor> &exec, QRSolver::GinkgoMatrix *A_Q, QRSolver::GinkgoMatrix *R) final override;
  ~HipQRSolver();

private:
  // Handles for HIP
#ifdef __HIPCC__
hipsolverDnHandle_t solverHandle;
hipsolverStatus_t   hipsolverStatus;
hipError_t          hipErrorCode = hipSuccess;
#endif
};

#endif