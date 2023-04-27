#ifndef PRECICE_NO_MPI

#pragma once

#include "testing/TestContext.hpp"

using namespace precice;
using precice::testing::TestContext;

/// tests the vector gradient API functions writeBlockVectorGradient data and writeVectorGradientData
void testVectorGradientFunctions(const TestContext &context);

#endif
