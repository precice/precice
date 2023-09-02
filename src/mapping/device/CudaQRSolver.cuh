#ifdef PRECICE_WITH_CUDA
#pragma once

#include <ginkgo/ginkgo.hpp>

/**
 * Computes the QR decomposition using CUDA
*/
void computeQRDecompositionCuda(const int deviceId, const std::shared_ptr<gko::Executor> &exec, gko::matrix::Dense<> *A_Q, gko::matrix::Dense<> *R);

#endif
