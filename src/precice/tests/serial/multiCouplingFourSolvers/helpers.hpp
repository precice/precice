#ifndef PRECICE_NO_MPI

#pragma once

#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

using namespace precice;
using precice::testing::TestContext;

struct MultiCouplingFourSolversFixture : testing::WhiteboxAccessor {
  std::string _pathToTests;

  MultiCouplingFourSolversFixture();

  void reset();
};

void multiCouplingFourSolvers(const std::string configFile, const TestContext &context);

#endif