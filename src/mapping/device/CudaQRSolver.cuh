#ifdef PRECICE_WITH_CUDA
#pragma once

#include "mapping/GinkgoDefinitions"

/**
 * Computes the QR decomposition using CUDA
*/
void computeQRDecompositionCuda(const int deviceId, const std::shared_ptr<gko::Executor> &exec, GinkgoMatrix *A_Q, GinkgoVector *R);

#endif
