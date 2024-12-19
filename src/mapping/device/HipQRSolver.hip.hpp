#ifdef PRECICE_WITH_HIP
#pragma once

#include "mapping/GinkgoDefinitions.hpp"

/**
 * Computes the QR decomposition using Hip
*/
void computeQRDecompositionHip(const std::shared_ptr<gko::Executor> &exec, GinkgoMatrix *A_Q, GinkgoVector *R);

void solvewithQRDecompositionHip(const std::shared_ptr<gko::Executor> &exec, GinkgoMatrix *U, GinkgoVector *x, GinkgoVector *rhs, GinkgoMatrix *matQ, GinkgoVector *in_vec);
#endif
