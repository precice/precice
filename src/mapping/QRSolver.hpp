#pragma once

#include <ginkgo/ginkgo.hpp>

/**
 * Base class for the two QR decomposition implementations, which can run on the device.
 */
class QRSolver {
public:
  using GinkgoMatrix = gko::matrix::Dense<>;

  virtual void computeQR(const std::shared_ptr<gko::Executor> &exec, GinkgoMatrix *A_Q, GinkgoMatrix *R) = 0;

protected:
  // Important variables which track the state of the solver routines
  double *dTau           = nullptr;
  void *  dWork          = nullptr;
  void *  hWork          = nullptr;
  int *   devInfo        = nullptr;
  int     backupDeviceId = 0;
};
