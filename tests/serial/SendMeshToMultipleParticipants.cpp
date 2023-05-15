#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/Participant.hpp>
#include <vector>

/**
 * @brief Tests sending one mesh to multiple participants
 */
BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_CASE(SendMeshToMultipleParticipants)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank), "SolverThree"_on(1_rank));

  Eigen::Vector2d vertex{0.0, 0.0};
  double          value = 1.0;

  std::string meshName;
  if (context.isNamed("SolverOne")) {
    meshName = "MeshA";
  } else if (context.isNamed("SolverTwo")) {
    meshName = "MeshB";
  } else if (context.isNamed("SolverThree")) {
    meshName = "MeshC";
  }

  precice::Participant interface(context.name, context.config(), 0, 1);

  const precice::VertexID vertexID = interface.setMeshVertex(meshName, vertex);
  auto                    dataName = "Data";
  interface.initialize();
  double maxDt = interface.getMaxTimeStepSize();

  if (context.isNamed("SolverOne")) {
    interface.writeData(meshName, dataName, {&vertexID, 1}, {&value, 1});
  } else {
    double valueReceived = -1.0;
    interface.readData(meshName, dataName, {&vertexID, 1}, maxDt, {&valueReceived, 1});
    BOOST_TEST(valueReceived == value);
  }

  interface.advance(maxDt);
  interface.finalize();
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial

#endif // PRECICE_NO_MPI
