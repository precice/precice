#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <boost/test/tools/detail/per_element_manip.hpp>
#include <precice/Participant.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(Compositional)
PRECICE_TEST_SETUP("M1SM"_on(1_rank), "M2SM"_on(1_rank), "M1"_on(1_rank), "M2"_on(1_rank), "Tendon"_on(1_rank));
BOOST_AUTO_TEST_CASE(FiveParticipants)
{
  PRECICE_TEST();

  precice::Participant participant(context.name, context.config(), context.rank, context.size);

  const std::vector<double> surfaceCoords{1, 0, 2, 0};
  const std::vector<double> neuralCoords{0, 0};

  std::vector<int> surface1VertexIDs(2);
  std::vector<int> surface2VertexIDs(2);

  std::vector<int> activationVertexIDs(1);
  std::vector<int> stretchVertexIDs(1);

  double timestepSize = 1.0;

  if (context.isNamed("M1SM")) {

    participant.setMeshVertices("Surface_M1SM_Mesh", surfaceCoords, surface1VertexIDs);
    participant.setMeshVertices("Activation_M1SM_Mesh", neuralCoords, activationVertexIDs);
    participant.setMeshVertices("Stretch_M1SM_Mesh", neuralCoords, stretchVertexIDs);

  } else if (context.isNamed("Tendon")) {
    participant.setMeshVertices("SurfaceTendon_M1SM_Mesh", surfaceCoords, surface1VertexIDs);
    participant.setMeshVertices("SurfaceTendon_M2SM_Mesh", surfaceCoords, surface2VertexIDs);

  } else if (context.isNamed("M2SM")) {

    participant.setMeshVertices("Surface_M2SM_Mesh", surfaceCoords, surface2VertexIDs);
    participant.setMeshVertices("Activation_M2SM_Mesh", neuralCoords, activationVertexIDs);
    participant.setMeshVertices("Stretch_M2SM_Mesh", neuralCoords, stretchVertexIDs);

  } else if (context.isNamed("M1")) {

    participant.setMeshVertices("Stretch_M1_Mesh", neuralCoords, stretchVertexIDs);
    participant.setMeshVertices("Activation_M1_Mesh", neuralCoords, activationVertexIDs);

  } else {
    participant.setMeshVertices("Stretch_M2_Mesh", neuralCoords, stretchVertexIDs);
    participant.setMeshVertices("Activation_M2_Mesh", neuralCoords, activationVertexIDs);
  }

  participant.initialize();

  std::vector<double> displacements1{1.2, 1.4};
  std::vector<double> displacements2{2.2, 2.4};
  std::vector<double> tractions1{1.2, 1.4};
  std::vector<double> tractions2{2.2, 2.4};
  std::vector<double> activation1{1.0};
  std::vector<double> activation2{2.0};
  std::vector<double> stretch1{1.1};
  std::vector<double> stretch2{2.2};

  std::vector<double> receivedDisplacements1{0.0, 0.0};
  std::vector<double> receivedDisplacements2{0.0, 0.0};
  std::vector<double> receivedTractions1{0.0, 0.0};
  std::vector<double> receivedTractions2{0.0, 0.0};
  std::vector<double> receivedActivation1{0.0};
  std::vector<double> receivedActivation2{0.0};
  std::vector<double> receivedStretch1{0.0};
  std::vector<double> receivedCrossStretch1{0.0};
  std::vector<double> receivedStretch2{0.0};
  std::vector<double> receivedCrossStretch2{0.0};

  for (int timestep = 0; timestep < 2; ++timestep) {

    if (context.isNamed("M1SM")) {
      participant.writeData("Surface_M1SM_Mesh", "Displacement1", surface1VertexIDs, displacements1);
      participant.readData("Surface_M1SM_Mesh", "Traction1", surface1VertexIDs, timestepSize, tractions1);
      participant.readData("Activation_M1SM_Mesh", "Activation1", activationVertexIDs, timestepSize, receivedActivation1);
      participant.writeData("Stretch_M1SM_Mesh", "Stretch1", stretchVertexIDs, stretch1);

    } else if (context.isNamed("Tendon")) {
      participant.readData("SurfaceTendon_M1SM_Mesh", "Displacement1", surface1VertexIDs, timestepSize, receivedDisplacements1);
      participant.readData("SurfaceTendon_M2SM_Mesh", "Traction2", surface2VertexIDs, timestepSize, receivedTractions2);
      participant.writeData("SurfaceTendon_M2SM_Mesh", "Displacement2", surface2VertexIDs, displacements2);
      participant.writeData("SurfaceTendon_M1SM_Mesh", "Traction1", surface1VertexIDs, tractions1);
    } else if (context.isNamed("M2SM")) {

      participant.readData("Surface_M2SM_Mesh", "Displacement2", surface2VertexIDs, timestepSize, receivedDisplacements2);
      participant.writeData("Surface_M2SM_Mesh", "Traction2", surface2VertexIDs, tractions2);
      participant.writeData("Stretch_M2SM_Mesh", "Stretch2", stretchVertexIDs, stretch2);

      participant.readData("Activation_M2SM_Mesh", "Activation2", activationVertexIDs, timestepSize, receivedActivation2);

    } else if (context.isNamed("M2")) {
      participant.readData("Stretch_M2_Mesh", "Stretch1", stretchVertexIDs, timestepSize, receivedCrossStretch2);
      participant.readData("Stretch_M2_Mesh", "Stretch2", stretchVertexIDs, timestepSize, receivedStretch2);
      participant.writeData("Activation_M2_Mesh", "Activation2", stretchVertexIDs, activation2);

    } else {

      BOOST_TEST(context.isNamed("M1"));

      participant.writeData("Activation_M1_Mesh", "Activation1", activationVertexIDs, activation1);
      participant.readData("Stretch_M1_Mesh", "Stretch1", stretchVertexIDs, timestepSize, receivedStretch1);
      participant.readData("Stretch_M1_Mesh", "Stretch2", stretchVertexIDs, timestepSize, receivedCrossStretch1);
    }

    if (participant.requiresWritingCheckpoint()) {
    }
    participant.advance(timestepSize);
    if (participant.requiresReadingCheckpoint()) {
    }
  }

  // Test read and write
  if (context.isNamed("M1SM")) {

    BOOST_TEST(receivedActivation1 == activation1, boost::test_tools::per_element());

  } else if (context.isNamed("Tendon")) {
    BOOST_TEST(receivedDisplacements1 == displacements1, boost::test_tools::per_element());

    BOOST_TEST(receivedTractions2 == tractions2, boost::test_tools::per_element());

  } else if (context.isNamed("M2SM")) {

    BOOST_TEST(receivedDisplacements2 == displacements2, boost::test_tools::per_element());

  } else if (context.isNamed("M1")) {

    BOOST_TEST(receivedStretch1 == stretch1, boost::test_tools::per_element());
    BOOST_TEST(receivedCrossStretch1 == stretch2, boost::test_tools::per_element());
  } else {
    BOOST_TEST(receivedCrossStretch2 == stretch1, boost::test_tools::per_element());
  }

  participant.finalize();
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Compositional

#endif // PRECICE_NO_MPI
