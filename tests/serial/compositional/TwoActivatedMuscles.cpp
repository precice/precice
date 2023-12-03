#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <boost/test/tools/detail/per_element_manip.hpp>
#include <precice/Participant.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(Compositional)
BOOST_AUTO_TEST_CASE(TwoActivatedMuscles)
{
  PRECICE_TEST("M1SM"_on(1_rank), "M2SM"_on(1_rank), "M1"_on(1_rank), "M2"_on(1_rank));

  precice::Participant participant(context.name, context.config(), context.rank, context.size);

  const std::vector<double> surfaceCoords{1, 0, 2, 0};
  const std::vector<double> neuralCoords{0, 0};

  std::vector<int> surfaceVertexIDs(2);
  std::vector<int> activationVertexIDs(1);
  std::vector<int> stretchVertexIDs(1);

  double timestepSize = 1.0;

  if (context.isNamed("M1SM")) {

    participant.setMeshVertices("Surface_M1SM_Mesh", surfaceCoords, surfaceVertexIDs);
    participant.setMeshVertices("Activation_M1SM_Mesh", neuralCoords, activationVertexIDs);
    participant.setMeshVertices("Stretch_M1SM_Mesh", neuralCoords, stretchVertexIDs);

  } else if (context.isNamed("M2SM")) {

    participant.setMeshVertices("Surface_M2SM_Mesh", surfaceCoords, surfaceVertexIDs);
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

  std::vector<double> tractions{1.2, 3.4};
  std::vector<double> displacements{4.2, 1.4};
  std::vector<double> activation1{1.0};
  std::vector<double> activation2{2.0};
  std::vector<double> stretch1{1.1};
  std::vector<double> stretch2{2.2};

  std::vector<double> receivedDisplacements{0.0, 0.0};
  std::vector<double> receivedActivation1{0.0};
  std::vector<double> receivedActivation2{0.0};
  std::vector<double> receivedStretch1{0.0};
  std::vector<double> receivedCrossStretch1{0.0};
  std::vector<double> receivedStretch2{0.0};

  for (int timestep = 0; timestep < 2; ++timestep) {

    if (context.isNamed("M1SM")) {
      participant.writeData("Surface_M1SM_Mesh", "Displacement", surfaceVertexIDs, displacements);
      participant.readData("Activation_M1SM_Mesh", "Activation1", activationVertexIDs, timestepSize, receivedActivation1);
      participant.writeData("Stretch_M1SM_Mesh", "stretch1", stretchVertexIDs, stretch1);

    } else if (context.isNamed("M2SM")) {

      participant.readData("Surface_M2SM_Mesh", "Displacement", surfaceVertexIDs, timestepSize, receivedDisplacements);
      participant.readData("Activation_M2SM_Mesh", "Activation2", activationVertexIDs, timestepSize, receivedActivation2);
      participant.writeData("Stretch_M2SM_Mesh", "stretch2", stretchVertexIDs, stretch2);

    } else if (context.isNamed("M1")) {

      participant.writeData("Activation_M1_Mesh", "Activation1", activationVertexIDs, activation1);
      participant.readData("Stretch_M1_Mesh", "stretch1", stretchVertexIDs, timestepSize, receivedStretch1);
      participant.readData("Stretch_M1_Mesh", "stretch2", stretchVertexIDs, timestepSize, receivedCrossStretch1);

    } else {

      BOOST_TEST(context.isNamed("M2"));

      participant.writeData("Activation_M2_Mesh", "Activation2", activationVertexIDs, activation2);
      participant.readData("Stretch_M2_Mesh", "stretch2", stretchVertexIDs, timestepSize, receivedStretch2);
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

  } else if (context.isNamed("M2SM")) {

    BOOST_TEST(receivedDisplacements == displacements, boost::test_tools::per_element());
    BOOST_TEST(receivedActivation2 == activation2, boost::test_tools::per_element());

  } else if (context.isNamed("M1")) {

    BOOST_TEST(receivedStretch1 == stretch1, boost::test_tools::per_element());
    BOOST_TEST(receivedCrossStretch1 == stretch2, boost::test_tools::per_element());

  } else {

    BOOST_TEST(receivedStretch2 == stretch2, boost::test_tools::per_element());
  }

  participant.finalize();
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Compositional

#endif // PRECICE_NO_MPI
