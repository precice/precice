#ifdef PRECICE_WITH_CUDA
#pragma once

#include "mapping/GinkgoDefinitions.hpp"

/**
 * Computes the QR decomposition using CUDA
*/
void computeQRDecompositionCuda(const std::shared_ptr<gko::Executor> &exec, GinkgoMatrix *A_Q, GinkgoVector *R);

void solvewithQRDecompositionCuda(const std::shared_ptr<gko::Executor> &exec, gko::matrix::Dense<> *U, gko::matrix::Dense<> *x, gko::matrix::Dense<> *rhs, gko::matrix::Dense<> *matQ,  gko::matrix::Dense<> *in_vec);

#endif
