#ifndef PRECICE_NO_MPI

#pragma once

#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

struct WhiteboxTestFixture : precice::testing::WhiteboxAccessor {
  WhiteboxTestFixture() {}
};

#endif