#ifndef PRECICE_NO_MPI

#pragma once

#include "testing/TestContext.hpp"

using namespace precice;
using precice::testing::TestContext;

void runTestQN(std::string const &config, TestContext const &context);

void runTestQNEmptyPartition(std::string const &config, TestContext const &context);

void runTestQNBoundedValue(std::string const &config, TestContext const &context);

#endif
