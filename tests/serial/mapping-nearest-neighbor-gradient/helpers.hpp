#ifndef PRECICE_NO_MPI

#pragma once

#include "testing/TestContext.hpp"

using namespace precice;
using precice::testing::TestContext;

/// tests the vector gradient API functions writeBlockVectorGradient data and writeVectorGradientData with the additional columnsFirst argument
void testVectorGradientFunctions(const TestContext &context, const bool writeBlockWise, const bool columnsFirst);

#endif