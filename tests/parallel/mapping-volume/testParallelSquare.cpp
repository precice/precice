#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Parallel)
BOOST_AUTO_TEST_SUITE(MappingVolume)
BOOST_AUTO_TEST_CASE(testParallelSquare)
{
  PRECICE_TEST("SolverOne"_on(2_ranks), "SolverTwo"_on(2_ranks));

  using precice::testing::equals;
  using precice::VertexID;

  // Implement your test here.
  BOOST_TEST(true);
  precice::SolverInterface interface(context.name, context.config(), context.rank, context.size);

  std::vector<VertexID> vertexIDs;

  if (context.isNamed("SolverOne")) {
    auto meshID = interface.getMeshID("MeshOne");
    auto dataID = interface.getDataID("DataOne", meshID);

    std::vector<double> coords;

    // Create a square with top left corner (rank 0) or bottom right. Diagonal "y = x" is shared.
    if (context.rank == 0) {
      coords = {0.0, 0.0, 
                1.0, 1.0, 
                0.0, 1.0};
    }
    else {
      coords = {0.0, 0.0, 
                1.0, 1.0, 
                1.0, 0.0};
    }

    interface.setMeshVertices(meshID, 3, coords.data(), vertexIDs.data());
    interface.setMeshTriangleWithEdges(meshID, vertexIDs[0], vertexIDs[1], vertexIDs[2]);
  } else {
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Parallel
BOOST_AUTO_TEST_SUITE_END() // MappingVolume

#endif // PRECICE_NO_MPI
