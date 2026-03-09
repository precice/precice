#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(Configuration)

/// Regression test for #1836: receive data declared but not exchanged should raise
/// a clear error instead of an assertion when mapping or reading.
PRECICE_TEST_SETUP("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank))
BOOST_AUTO_TEST_CASE(ReceiveDataNotExchanged)
{
  PRECICE_TEST();
  precice::Participant interface(context.name, context.config(), context.rank, context.size);

  if (context.isNamed("SolverOne")) {
    auto   meshName = "MeshOne";
    double coords[] = {0.1, 1.2};
    auto   vertexid = interface.setMeshVertex(meshName, coords);

    double dataOne[] = {3.4};
    interface.writeData(meshName, "DataOne", {&vertexid, 1}, dataOne);
  } else {
    auto   meshName = "MeshTwo";
    double coords[] = {0.12, 1.21};
    auto   vertexid = interface.setMeshVertex(meshName, coords);

    double dataTwo[] = {7.8};
    interface.writeData(meshName, "DataTwo", {&vertexid, 1}, dataTwo);
  }

  if (context.isNamed("SolverOne")) {
    try {
      interface.initialize();
      double            dt       = interface.getMaxTimeStepSize();
      auto              meshName = "MeshOne";
      double            dataTwo[1];
      precice::VertexID vid = 0;
      interface.readData(meshName, "DataTwo", {&vid, 1}, 0.0, dataTwo);
      double dataOne[] = {3.4};
      interface.writeData(meshName, "DataOne", {&vid, 1}, dataOne);
      interface.advance(dt);
    } catch (const ::precice::Error &) {
      // SolverTwo failed during init (receive data not exchanged); partner disconnected
    }
  } else {
    BOOST_CHECK_EXCEPTION(
        interface.initialize(),
        ::precice::Error,
        precice::testing::errorContains("contain any data samples"));
  }
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

#endif
