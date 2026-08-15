#ifndef PRECICE_NO_MPI

#include <array>
#include <boost/test/unit_test_log.hpp>
#include <precice/precice.hpp>

#include "helpers.hpp"
#include "testing/Testing.hpp"

/// Test to run simple "do nothing" coupling between two solvers.
void runExporter(std::string const &configurationFileName, precice::testing::TestContext const &context)
{
  BOOST_TEST_MESSAGE("Config: " << configurationFileName);
  precice::Participant p(context.name, configurationFileName, 0, 1);

  std::array<double, 6> pos{0, 0, 0, 1, 1, 1};
  std::array<int, 2>    vids;

  auto meshName = context.isNamed("SA") ? "MA" : "MB";
  BOOST_REQUIRE(p.getMeshDimensions(meshName) == 3);
  p.setMeshVertices(meshName, pos, vids);

  p.initialize();
  while (p.isCouplingOngoing()) {
    p.advance(p.getMaxTimeStepSize());
  }
}

#endif
