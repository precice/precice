#ifndef PRECICE_NO_MPI

#pragma once

#include "testing/TestContext.hpp"

using namespace precice;
using precice::testing::TestContext;

void runTestAccessReceivedMesh(const TestContext &       context,
                               const std::vector<double> boundingBoxSecondary,
                               const std::vector<double> writeDataSecondary,
                               const std::vector<double> expectedPositionSecondary,
                               const std::vector<double> expectedReadDataSecondary,
                               const int                 startIndex);

#endif