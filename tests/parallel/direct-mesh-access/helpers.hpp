#ifndef PRECICE_NO_MPI

#pragma once

#include "testing/TestContext.hpp"

using namespace precice;
using precice::testing::TestContext;

void runTestAccessReceivedMesh(const TestContext &context, const std::vector<double> boundingBoxSlave,
                               const std::vector<double> writeDataSlave,
                               const std::vector<double> expectedPositionSlave,
                               const std::vector<double> expectedReadDataSlave, const int startIndex);

#endif