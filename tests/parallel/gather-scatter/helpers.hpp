#ifndef PRECICE_NO_MPI

#pragma once

#include "testing/TestContext.hpp"

using namespace precice;
using precice::testing::TestContext;

void runTestEnforceGatherScatter(std::vector<double> masterPartition, const TestContext &context);

#endif