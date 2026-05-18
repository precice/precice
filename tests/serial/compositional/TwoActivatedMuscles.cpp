#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <boost/test/tools/detail/per_element_manip.hpp>
#include <precice/Participant.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(Compositional)
PRECICE_TEST_SETUP("M1SM"_on(1_rank), "M2SM"_on(1_rank), "M1"_on(1_rank), "M2"_on(1_rank))
BOOST_AUTO_TEST_CASE(TwoActivatedMuscles)
{
  PRECICE_TEST();

  precice::Participant participant(context.name, context.config(), context.rank, context.size);

  const std::vector<double> surfaceCoords{1, 0, 2, 0};
  const std::vector<double> neuralCoords{0, 0};

  std::vector<int> surfaceVertexIDs(2);
  std::vector<int> activationVertexIDs(1);
  std::vector<int> stretchVertexIDs(1);

  const double timestepSize = 1.0;

  if (context.isNamed("M1SM")) {
    participant.setMeshVertices("Surface-M1SM-Mesh", surfaceCoords, surfaceVertexIDs);
    participant.setMeshVertices("Activation-M1SM-Mesh", neuralCoords, activationVertexIDs);
    participant.setMeshVertices("Stretch-M1SM-Mesh", neuralCoords, stretchVertexIDs);
  }
  if (context.isNamed("M2SM")) {
    participant.setMeshVertices("Surface-M2SM-Mesh", surfaceCoords, surfaceVertexIDs);
    participant.setMeshVertices("Activation-M2SM-Mesh", neuralCoords, activationVertexIDs);
    participant.setMeshVertices("Stretch-M2SM-Mesh", neuralCoords, stretchVertexIDs);
  }
  if (context.isNamed("M1")) {
    participant.setMeshVertices("Stretch-M1-Mesh", neuralCoords, stretchVertexIDs);
    participant.setMeshVertices("Activation-M1-Mesh", neuralCoords, activationVertexIDs);
  }
  if (context.isNamed("M2")) {
    participant.setMeshVertices("Stretch-M2-Mesh", neuralCoords, stretchVertexIDs);
    participant.setMeshVertices("Activation-M2-Mesh", neuralCoords, activationVertexIDs);
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

  const bool isImplicit = context.isNamed("M1SM") || context.isNamed("M2SM");
  for (int timestep = 0; timestep < 2; ++timestep) {
    BOOST_REQUIRE(participant.isCouplingOngoing());

    if (isImplicit) {
      BOOST_TEST(participant.requiresWritingCheckpoint());
    } else {
      BOOST_TEST(!participant.requiresWritingCheckpoint());
    }

    if (context.isNamed("M1SM")) {
      participant.writeData("Surface-M1SM-Mesh", "Displacement", surfaceVertexIDs, displacements);
      participant.readData("Activation-M1SM-Mesh", "Activation-M1", activationVertexIDs, timestepSize, receivedActivation1);
      participant.writeData("Stretch-M1SM-Mesh", "Stretch-M1SM", stretchVertexIDs, stretch1);
    }
    if (context.isNamed("M2SM")) {
      participant.readData("Surface-M2SM-Mesh", "Displacement", surfaceVertexIDs, timestepSize, receivedDisplacements);
      participant.readData("Activation-M2SM-Mesh", "Activation-M2", activationVertexIDs, timestepSize, receivedActivation2);
      participant.writeData("Stretch-M2SM-Mesh", "Stretch-M2SM", stretchVertexIDs, stretch2);
    }
    if (context.isNamed("M1")) {
      participant.writeData("Activation-M1-Mesh", "Activation-M1", activationVertexIDs, activation1);
      participant.readData("Stretch-M1-Mesh", "Stretch-M1SM", stretchVertexIDs, timestepSize, receivedStretch1);
      participant.readData("Stretch-M1-Mesh", "Stretch-M2SM", stretchVertexIDs, timestepSize, receivedCrossStretch1);
    }
    if (context.isNamed("M2")) {
      participant.writeData("Activation-M2-Mesh", "Activation-M2", activationVertexIDs, activation2);
      participant.readData("Stretch-M2-Mesh", "Stretch-M2SM", stretchVertexIDs, timestepSize, receivedStretch2);
    }

    participant.advance(timestepSize);

    // Always converges
    BOOST_TEST(!participant.requiresReadingCheckpoint());
  }

  BOOST_REQUIRE(!participant.isCouplingOngoing());

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
