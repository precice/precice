#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <boost/test/tools/detail/per_element_manip.hpp>
#include <precice/Participant.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(Compositional)
BOOST_AUTO_TEST_CASE(OneActivatedMuscle)
{
  PRECICE_TEST("M1SM"_on(1_rank), "M2SM"_on(1_rank), "Tendon"_on(1_rank), "M1"_on(1_rank));

  precice::Participant participant(context.name, context.config(), context.rank, context.size);

  const std::vector<double> surfaceCoords{1, 0, 2, 0};
  const std::vector<double> neuralCoords{0, 0};

  std::vector<int> surface1VertexIDs(2);
  std::vector<int> surface2VertexIDs(2);
  std::vector<int> activationVertexIDs(1);
  std::vector<int> stretchVertexIDs(1);

  const double timestepSize = 1.0;

  if (context.isNamed("M1")) {
    participant.setMeshVertices("Activation_M1_Mesh", neuralCoords, activationVertexIDs);
    participant.setMeshVertices("Stretch_M1_Mesh", neuralCoords, stretchVertexIDs);
  }
  if (context.isNamed("M1SM")) {
    participant.setMeshVertices("Surface_M1SM_Mesh", surfaceCoords, surface1VertexIDs);
    participant.setMeshVertices("Activation_M1SM_Mesh", neuralCoords, activationVertexIDs);
    participant.setMeshVertices("Stretch_M1SM_Mesh", neuralCoords, stretchVertexIDs);
  }
  if (context.isNamed("M2SM")) {
    participant.setMeshVertices("Surface_M2SM_Mesh", surfaceCoords, surface2VertexIDs);
    participant.setMeshVertices("Stretch_M2SM_Mesh", neuralCoords, stretchVertexIDs);
  }
  if (context.isNamed("Tendon")) {
    participant.setMeshVertices("SurfaceTendon_M1SM_Mesh", surfaceCoords, surface1VertexIDs);
    participant.setMeshVertices("SurfaceTendon_M2SM_Mesh", surfaceCoords, surface2VertexIDs);
  }

  std::vector<double> activation1{1.0};
  std::vector<double> stretch1{1.1};
  std::vector<double> stretch2{2.2};
  std::vector<double> tractions1{1.2, 3.4};
  std::vector<double> displacements1{4.2, 1.4};
  std::vector<double> tractions2{1.2, 3.7};
  std::vector<double> displacements2{4.1, 1.4};

  std::vector<double> receivedDisplacements{0.0, 0.0};
  std::vector<double> receivedActivation1{0.0};
  std::vector<double> receivedStretch1{0.0};
  std::vector<double> receivedStretch2{0.0};

  participant.initialize();

  bool isImplicit = context.isNamed("Tendon") || context.isNamed("M1SM") || context.isNamed("M2SM");

  for (int timestep = 1; timestep <= 2; ++timestep) {
    BOOST_REQUIRE(participant.isCouplingOngoing()); // Tendon M2SM

    if (isImplicit) {
      BOOST_TEST(participant.requiresWritingCheckpoint());
    } else {
      BOOST_TEST(!participant.requiresWritingCheckpoint());
    }

    if (context.isNamed("M1")) {
      participant.writeData("Activation_M1_Mesh", "Activation1", activationVertexIDs, activation1);
      participant.readData("Stretch_M1_Mesh", "stretch1", stretchVertexIDs, timestepSize, receivedStretch1);
      participant.readData("Stretch_M1_Mesh", "stretch2", stretchVertexIDs, timestepSize, receivedStretch2);
    }
    if (context.isNamed("M1SM")) {
      participant.readData("Activation_M1SM_Mesh", "Activation1", activationVertexIDs, timestepSize, receivedActivation1);
      participant.writeData("Stretch_M1SM_Mesh", "stretch1", stretchVertexIDs, stretch1);
      participant.writeData("Surface_M1SM_Mesh", "Displacement1", surface1VertexIDs, displacements1);
    }
    if (context.isNamed("Tendon")) {
      participant.readData("SurfaceTendon_M1SM_Mesh", "Displacement1", surface1VertexIDs, timestepSize, receivedDisplacements);
    }
    if (context.isNamed("M2SM")) {
      participant.writeData("Stretch_M2SM_Mesh", "stretch2", stretchVertexIDs, stretch2);
    }

    participant.advance(timestepSize);

    // Schemes always converge
    BOOST_TEST(!participant.requiresReadingCheckpoint());
  }

  BOOST_REQUIRE(!participant.isCouplingOngoing());

  // Test read and write
  if (context.isNamed("Tendon")) {
    BOOST_TEST(receivedDisplacements == displacements1, boost::test_tools::per_element());
  }
  if (context.isNamed("M1SM")) {
    BOOST_TEST(receivedActivation1 == activation1, boost::test_tools::per_element());
  }
  if (context.isNamed("M1")) {
    BOOST_TEST(receivedStretch1 == stretch1, boost::test_tools::per_element());
    BOOST_TEST(receivedStretch2 == stretch2, boost::test_tools::per_element());
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Compositional

#endif // PRECICE_NO_MPI
