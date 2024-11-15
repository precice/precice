#ifndef PRECICE_NO_MPI

#pragma once

#include "testing/TestContext.hpp"

using namespace precice;
using precice::testing::TestContext;

void subcyclingWithNSteps(TestContext const &context, int nSubsteps, bool useAdvancedDtStrategy = false);

#endif
