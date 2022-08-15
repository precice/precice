#ifndef PRECICE_NO_MPI

#pragma once

#include "testing/TestContext.hpp"

using namespace precice;
using precice::testing::TestContext;

/// Tests stationary mapping with solver provided meshes.
void runTestStationaryMappingWithSolverMesh(std::string const &config, int dim, TestContext const &context);

#endif
