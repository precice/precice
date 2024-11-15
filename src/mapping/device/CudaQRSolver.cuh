#ifdef PRECICE_WITH_CUDA
#pragma once

#include "mapping/GinkgoDefinitions.hpp"

/**
 * Computes the QR decomposition using CUDA
*/
void computeQRDecompositionCuda(const std::shared_ptr<gko::Executor> &exec, GinkgoMatrix *A_Q, GinkgoVector *R);

#endif
