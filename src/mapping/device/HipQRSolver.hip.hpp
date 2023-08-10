#ifdef PRECICE_WITH_HIP
#pragma once

#include <ginkgo/ginkgo.hpp>

/**
 * Computes the QR decomposition using Hip
*/
void computeQRDecompositionHip(const int deviceId, const std::shared_ptr<gko::Executor> &exec, gko::matrix::Dense<> *A_Q, gko::matrix::Dense<> *R);

#endif
