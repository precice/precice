#ifndef PRECICE_NO_MPI

#pragma once

#include "testing/TestContext.hpp"

using namespace precice;
using precice::testing::TestContext;

void testMappingVolumeOneTriangle(const std::string configFile, const TestContext &context, bool read);
void testMappingVolumeOneTriangleConservative(const std::string configFile, const TestContext &context);
void testMappingVolumeOneTetra(const std::string configFile, const TestContext &context, bool read);
void testMappingVolumeOneTetraConservative(const std::string configFile, const TestContext &context);

#endif
