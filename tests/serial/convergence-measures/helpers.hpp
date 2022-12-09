#ifndef PRECICE_NO_MPI

#pragma once

#include "testing/TestContext.hpp"

using namespace precice;
using precice::testing::TestContext;

/**
 * @brief method to test whether certain convergence measures give the correct number of iterations
 */
void testConvergenceMeasures(const std::string configFile, TestContext const &context, std::vector<int> &expectedIterations);

#endif
