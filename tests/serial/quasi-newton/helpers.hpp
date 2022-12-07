#ifndef PRECICE_NO_MPI

#pragma once

#include "testing/TestContext.hpp"

using namespace precice;
using precice::testing::TestContext;

namespace tests::serial {
void runTestQN(std::string const &config, TestContext const &context);
} // namespace tests::serial

#endif
