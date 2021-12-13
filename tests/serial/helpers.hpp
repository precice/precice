#ifndef PRECICE_NO_MPI

#pragma once

#include "testing/TestContext.hpp"

using namespace precice;
using precice::testing::TestContext;

void multiCouplingFourSolvers(const std::string configFile, const TestContext &context);

#endif