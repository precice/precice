#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

// Test representing the full explicit lifecycle of a Participant
BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Parallel)
BOOST_AUTO_TEST_SUITE(Lifecycle)
BOOST_AUTO_TEST_CASE(Full)
{
  PRECICE_TEST("SolverOne"_on(2_ranks), "SolverTwo"_on(2_ranks));

  precice::Participant interface(context.name, context.config(), context.rank, context.size);
  constexpr double     y{0};
  constexpr double     z{0};
  constexpr double     x1{1};
  constexpr double     dx{1};

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
  interface.finalize();
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Parallel
BOOST_AUTO_TEST_SUITE_END() // Lifecycle

#endif // PRECICE_NO_MPI
