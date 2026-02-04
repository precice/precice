#include <boost/test/tools/interface.hpp>
#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(Implicit)
PRECICE_TEST_SETUP("PA"_on(1_rank), "PB"_on(1_rank))
BOOST_AUTO_TEST_CASE(NoRead)
{
  PRECICE_TEST();
  /// Test simple coupled simulation with coupling iterations.

  std::string meshName  = context.isNamed("PA") ? "MA" : "MB";
  std::string readName  = context.isNamed("PA") ? "DB" : "DA";
  std::string writeName = context.isNamed("PA") ? "DA" : "DB";

  precice::Participant p(context.name, context.config(), context.rank, context.size);

  std::array<double, 3>            pos{0.0, 0.0, 0.0};
  std::array<precice::VertexID, 1> vids;
  p.setMeshVertices(meshName, pos, vids);
  p.initialize();

  std::array<double, 1> data{0.0};

  BOOST_REQUIRE(p.isCouplingOngoing());

  double dt = p.getMaxTimeStepSize();
  BOOST_TEST(p.requiresWritingCheckpoint());
  p.readData(meshName, readName, vids, dt, data);
  data.fill(1.0);
  p.writeData(meshName, writeName, vids, data);
  p.advance(dt);
  // reading is required
  BOOST_CHECK_EXCEPTION(p.advance(dt),
                        ::precice::Error,
                        ::precice::testing::errorContains("requiresReadingCheckpoint()"));
}

BOOST_AUTO_TEST_SUITE_END() // Implicit
BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial

#endif // PRECICE_NO_MPI
