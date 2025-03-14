#ifndef PRECICE_NO_MPI

#pragma once

#include "testing/TestContext.hpp"

using namespace precice;
using precice::testing::TestContext;

std::vector<double> readDoublesFromTXTFile(const std::string &filename, int skip = 0);

#endif
