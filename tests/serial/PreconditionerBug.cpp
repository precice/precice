#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
#include <vector>

/**
 * @brief Test to reproduce the problem of issue 383, https://github.com/precice/precice/issues/383
 */
BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_CASE(PreconditionerBug)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  using Eigen::Vector2d;

  std::string              meshName = context.isNamed("SolverOne") ? "MeshOne" : "MeshTwo";
  precice::SolverInterface interface(context.name, context.config(), context.rank, context.size);

  Vector2d vertex{0.0, 0.0};

  precice::VertexID vertexID = interface.setMeshVertex(meshName, vertex.data());

  interface.initialize();
  int numberOfAdvanceCalls = 0;

  while (interface.isCouplingOngoing()) {
    if (interface.requiresWritingCheckpoint()) {
      // nothing
    }
    if (context.isNamed("SolverTwo")) {
      auto dataName = "DataOne";
      // to get convergence in first timestep (everything 0), but not in second timestep
      Vector2d value{0.0, 2.0 + numberOfAdvanceCalls * numberOfAdvanceCalls};
      interface.writeVectorData(meshName, dataName, vertexID, value.data());
    }
    interface.advance(1.0);

    if (interface.requiresReadingCheckpoint()) {
      // nothing
    }
    ++numberOfAdvanceCalls;
  }
  interface.finalize();
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial

#endif // PRECICE_NO_MPI
