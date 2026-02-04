#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Fundamental)
BOOST_AUTO_TEST_SUITE(JustInTime)
BOOST_AUTO_TEST_SUITE(Both)

PRECICE_TEST_SETUP("Provider"_on(1_rank), "Receiver"_on(1_rank))
BOOST_AUTO_TEST_CASE(NoneConfigured)
{
  PRECICE_TEST();

  // Set up Participant
  precice::Participant p(context.name, context.config(), 0, 1);

  std::vector<double> positions = {0.0, 0.0, 1.0, 0.0};
  std::vector<int>    ids(2, -1);
  std::vector<double> data(2);

  if (context.isNamed("Provider")) {
    p.setMeshVertices("M", positions, ids);
    p.initialize();
    p.advance(p.getMaxTimeStepSize());
    BOOST_REQUIRE(!p.isCouplingOngoing());
  } else {
    std::array<double, 4> boundingBox = {0.0, 1.0, 0.0, 1.0};
    p.setMeshAccessRegion("M", boundingBox);
    p.initialize();

    std::vector<double> locations = {0.25, 0.0, 0.5, 0.0, 0.75, 0.0};
    std::vector<double> readData(3, 0.0);
    std::vector<double> writeData(3, 3.1415);
    auto                dt = p.getMaxTimeStepSize();
    BOOST_CHECK_EXCEPTION(p.mapAndReadData("M", "D1", locations, dt, readData), precice::Error, ::precice::testing::errorContains("no matching just-in-time mapping"));
    BOOST_CHECK_EXCEPTION(p.writeAndMapData("M", "D2", locations, writeData), precice::Error, ::precice::testing::errorContains("no matching just-in-time mapping"));

    p.advance(dt);
    BOOST_REQUIRE(!p.isCouplingOngoing());
  }
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

#endif // PRECICE_NO_MPI
