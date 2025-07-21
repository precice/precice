#ifndef PRECICE_NO_MPI

#pragma once

#include "testing/TestContext.hpp"

using namespace precice;
using precice::testing::TestContext;

void runTestAccessReceivedMesh(const TestContext        &context,
                               const std::vector<double> boundingBoxSecondaryRank,
                               const std::vector<double> writeDataSecondaryRank,
                               const std::vector<double> expectedPositionSecondaryRank,
                               const std::vector<double> expectedReadDataSecondaryRank,
                               const size_t              startIndex);

void runTestMultipleBoundingBoxes2D(const TestContext &context);

#endif
