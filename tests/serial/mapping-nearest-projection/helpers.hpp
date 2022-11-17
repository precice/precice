#ifndef PRECICE_NO_MPI

#pragma once

#include "testing/TestContext.hpp"

using namespace precice;
using precice::testing::TestContext;

void testMappingNearestProjection(bool defineEdgesExplicitly, bool useBulkFunctions, const std::string configFile, const TestContext &context);

void testQuadMappingNearestProjection(bool defineEdgesExplicitly, bool useBulkFunctions, const std::string configFile, const TestContext &context);

void testQuadMappingNearestProjectionTallKite(bool defineEdgesExplicitly, bool useBulkFunctions, const std::string configFile, const TestContext &context);

void testQuadMappingNearestProjectionWideKite(bool defineEdgesExplicitly, bool useBulkFunctions, const std::string configFile, const TestContext &context);

inline void testMappingNearestProjection(bool defineEdgesExplicitly, const std::string configFile, const TestContext &context)
{
  testMappingNearestProjection(defineEdgesExplicitly, false, configFile, context);
}

inline void testQuadMappingNearestProjection(bool defineEdgesExplicitly, const std::string configFile, const TestContext &context)
{
  testQuadMappingNearestProjection(defineEdgesExplicitly, false, configFile, context);
}

inline void testQuadMappingNearestProjectionTallKite(bool defineEdgesExplicitly, const std::string configFile, const TestContext &context)
{
  testQuadMappingNearestProjectionTallKite(defineEdgesExplicitly, false, configFile, context);
}

inline void testQuadMappingNearestProjectionWideKite(bool defineEdgesExplicitly, const std::string configFile, const TestContext &context)
{
  testQuadMappingNearestProjectionWideKite(defineEdgesExplicitly, false, configFile, context);
}

#endif
