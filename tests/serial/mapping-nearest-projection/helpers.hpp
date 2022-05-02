#ifndef PRECICE_NO_MPI

#pragma once

#include "testing/TestContext.hpp"

using namespace precice;
using precice::testing::TestContext;

void testMappingNearestProjection(bool defineEdgesExplicitly, const std::string configFile, const TestContext &context);

void testQuadMappingNearestProjection(bool defineEdgesExplicitly, const std::string configFile,
                                      const TestContext &context);

void testQuadMappingNearestProjectionTallKite(bool defineEdgesExplicitly, const std::string configFile,
                                              const TestContext &context);

void testQuadMappingNearestProjectionWideKite(bool defineEdgesExplicitly, const std::string configFile,
                                              const TestContext &context);

#endif