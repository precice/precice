#include <boost/test/tools/interface.hpp>
#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(Implicit)
PRECICE_TEST_SETUP("PA"_on(1_rank), "PB"_on(1_rank))
BOOST_AUTO_TEST_CASE(Basic)
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

  while (p.isCouplingOngoing()) {
    double dt = p.getMaxTimeStepSize();

    if (p.requiresWritingCheckpoint()) {
      data.fill(0.0); // reset
    } else {
      p.readData(meshName, readName, vids, dt, data);
    }

    // PA forwards the data
    if (context.isNamed("PB")) {
      double y = (data.front() + dt) / (1.0 + dt);
      data.fill(y);
    }
    p.writeData(meshName, writeName, vids, data);
    p.advance(dt); // converges in 7 iterations

    if (p.requiresReadingCheckpoint()) {
    }
  }
}

BOOST_AUTO_TEST_SUITE_END() // Implicit
BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial

#endif // PRECICE_NO_MPI
