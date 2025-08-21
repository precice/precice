#include <boost/test/tools/interface.hpp>
#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(Implicit)
PRECICE_TEST_SETUP("PA"_on(1_rank), "PB"_on(1_rank))
BOOST_AUTO_TEST_CASE(NoInit)
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
  BOOST_CHECK_EXCEPTION(p.initialize(),
                        ::precice::Error,
                        ::precice::testing::errorContains("requiresInitialData()"));
}

BOOST_AUTO_TEST_SUITE_END() // Implicit
BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial

#endif // PRECICE_NO_MPI
