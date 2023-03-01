#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
#include <vector>

// Test representing the full explicit lifecycle of a SolverInterface
BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Parallel)
BOOST_AUTO_TEST_SUITE(Lifecycle)
BOOST_AUTO_TEST_CASE(Full)
{
  PRECICE_TEST("SolverOne"_on(2_ranks), "SolverTwo"_on(2_ranks));

  precice::SolverInterface interface(context.name, context.config(), context.rank, context.size);
  constexpr double         y{0};
  constexpr double         z{0};
  constexpr double         x1{1};
  constexpr double         dx{1};

  if (context.isNamed("SolverOne")) {
    auto   meshid   = "MeshOne";
    double coords[] = {x1 + dx * context.rank, y, z};
    auto   vertexid = interface.setMeshVertex(meshid, coords);

    auto   dataid = "DataOne"; //  meshid
    double data[] = {3.4, 4.5, 5.6};
    interface.writeVectorData(meshID, dataid, vertexid, data);
  } else {
    auto   meshid   = "MeshTwo";
    double coords[] = {x1 + dx * context.rank, y, z};
    auto   vertexid = interface.setMeshVertex(meshid, coords);

    auto dataid = "DataTwo"; //  meshid
    interface.writeScalarData(meshID, dataid, vertexid, 7.8);
  }
  interface.initialize();
  BOOST_TEST(interface.isCouplingOngoing());
  interface.finalize();
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Parallel
BOOST_AUTO_TEST_SUITE_END() // Lifecycle

#endif // PRECICE_NO_MPI
