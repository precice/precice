#ifdef PRECICE_WITH_CUDA
#pragma once

#include "mapping/GinkgoDefinitions.hpp"

/**
 * Computes the QR decomposition using CUDA
*/
void computeQRDecompositionCuda(const std::shared_ptr<gko::Executor> &exec, GinkgoMatrix *A_Q, GinkgoVector *R);

void solvewithQRDecompositionCuda(const std::shared_ptr<gko::Executor> &exec, GinkgoMatrix *U, GinkgoVector *x, GinkgoVector *rhs, GinkgoMatrix *matQ,  GinkgoVector *in_vec);

#endif
