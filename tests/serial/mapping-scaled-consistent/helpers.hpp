#ifndef PRECICE_NO_MPI

#pragma once

#include "testing/TestContext.hpp"

using namespace precice;
using precice::testing::TestContext;

void testQuadMappingScaledConsistent(const std::string configFile, const TestContext &context);
void testQuadMappingScaledConsistentVolumetric(const std::string configFile, const TestContext &context);
void testTetraScaledConsistentVolumetric(const std::string configFile, const TestContext &context);

#endif
