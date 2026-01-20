#ifndef PRECICE_NO_MPI

#pragma once

#include "testing/TestContext.hpp"

using namespace precice;
using precice::testing::TestContext;

void testRBFMapping(const std::string configFile, const TestContext &context);
void testRBFMappingVectorial(const std::string configFile, const TestContext &context, bool mappingIsConservative);

#endif
