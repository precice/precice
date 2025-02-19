#ifndef PRECICE_NO_MPI

#pragma once

#include "testing/TestContext.hpp"

using namespace precice;
using precice::testing::TestContext;

void runTestQN(bool includeSecondaryData, std::string const &config, TestContext const &context);

void runTestQNEmptyPartition(std::string const &config, TestContext const &context);

void runTestQNWithWaveforms(std::string const &config, TestContext const &context);

void runTestQNWithWaveformsReducedTimeGrid(std::string const &config, TestContext const &context);
#endif
