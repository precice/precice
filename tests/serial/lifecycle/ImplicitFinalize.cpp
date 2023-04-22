#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(Lifecycle)
// Test representing the full lifecycle of a SolverInterface
// Finalize is not called explicitly here.
// The destructor has to cleanup.
BOOST_AUTO_TEST_CASE(ImplicitFinalize)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));
  precice::SolverInterface interface(context.name, context.config(), context.rank, context.size);

  if (context.isNamed("SolverOne")) {
    auto   meshName = "MeshOne";
    double coords[] = {0.1, 1.2, 2.3};
    auto   vertexid = interface.setMeshVertex(meshName, coords);

    auto   dataName = "DataOne";
    double data[]   = {3.4, 4.5, 5.6};
    interface.writeVectorData(meshName, dataName, vertexid, data);
  } else {
    auto   meshName = "MeshTwo";
    double coords[] = {0.12, 1.21, 2.2};
    auto   vertexid = interface.setMeshVertex(meshName, coords);

    auto dataName = "DataTwo";
    interface.writeScalarData(meshName, dataName, vertexid, 7.8);
  }
  interface.initialize();
  BOOST_TEST(interface.isCouplingOngoing());
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Lifecycle

#endif // PRECICE_NO_MPI
