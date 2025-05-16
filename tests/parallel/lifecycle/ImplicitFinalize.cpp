#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

// Test representing the full lifecycle of a Participant
// Finalize is not called explicitly here.
// The destructor has to cleanup.
BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Parallel)
BOOST_AUTO_TEST_SUITE(Lifecycle)
PRECICE_TEST_SETUP("SolverOne"_on(2_ranks), "SolverTwo"_on(2_ranks))
BOOST_AUTO_TEST_CASE(ImplicitFinalize)
{
  PRECICE_TEST();

  precice::Participant interface(context.name, context.config(), context.rank, context.size);

  constexpr double y{0};
  constexpr double z{0};
  constexpr double x1{1};
  constexpr double dx{1};

  if (context.isNamed("SolverOne")) {
    auto   meshName = "MeshOne";
    double coords[] = {x1 + dx * context.rank, y, z};
    auto   vertexid = interface.setMeshVertex(meshName, coords);

    auto   dataName = "DataOne";
    double data[]   = {3.4, 4.5, 5.6};
    interface.writeData(meshName, dataName, {&vertexid, 1}, data);
  } else {
    auto   meshName = "MeshTwo";
    double coords[] = {x1 + dx * context.rank, y, z};
    auto   vertexid = interface.setMeshVertex(meshName, coords);

    auto   dataName = "DataTwo";
    double data[]   = {7.8};
    interface.writeData(meshName, dataName, {&vertexid, 1}, data);
  }
  interface.initialize();
  BOOST_TEST(interface.isCouplingOngoing());
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Parallel
BOOST_AUTO_TEST_SUITE_END() // Lifecycle

#endif // PRECICE_NO_MPI
