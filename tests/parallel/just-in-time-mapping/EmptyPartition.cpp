#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Parallel)
BOOST_AUTO_TEST_SUITE(JustInTimeMapping)
PRECICE_TEST_SETUP("SolverOne"_on(2_ranks), "SolverTwo"_on(2_ranks))
BOOST_AUTO_TEST_CASE(EmptyPartition)
{
  PRECICE_TEST();
  // Test case for a just-in-time mapping, where the access region is empty on
  // a particular rank. Here, SolverOne does not receive any vertices due to the
  // defined bounding box on rank 1. Each solver runs on two ranks. SolverTwo defines
  // 5(2 and 3) vertices which need to be repartitioned on SolverOne according to the
  // defined boundingBoxes (resulting in 3  vertices on one rank and 2 completely filtered
  // vertices). Filtered vertices are filled with zero data values
  const std::vector<double> boundingBoxSecondaryRank      = std::vector<double>{10.0, 12.0, 12.0, 17};
  const std::vector<double> expectedPositionSecondaryRank = std::vector<double>{};
  const std::vector<double> expectedReadDataSecondaryRank = std::vector<double>({3., 0., 0.});

  if (context.isNamed("SolverOne")) {
    // Defines the bounding box and writes data to the received mesh
    precice::Participant interface(context.name, context.config(), context.rank, context.size);
    auto                 otherMeshName   = "MeshTwo";
    auto                 velocityName    = "Velocities";
    auto                 uncertaintyName = "Uncertainty";
    const int            dim             = interface.getMeshDimensions(otherMeshName);

    std::vector<double> boundingBox = context.isPrimary() ? std::vector<double>({0.0, 1.0, 0.0, 3.5}) : boundingBoxSecondaryRank;
    // Set bounding box
    interface.setMeshAccessRegion(otherMeshName, boundingBox);
    // Initialize the Participant
    interface.initialize();
    double dt = interface.getMaxTimeStepSize();

    // Allocate memory
    std::vector<double> positions = context.isPrimary() ? std::vector<double>({0.5, 2.0, 1.0, 3.5}) : std::vector<double>({11.0, 13.0, 12.0, 12.0});

    // Create some unique writeData in order to check it in the other participant
    std::vector<double> writeData(2, 5);

    while (interface.isCouplingOngoing()) {
      // Write data
      if (context.isPrimary()) {
        // interface.writeData(otherMeshName, velocityName, ids, writeData);
        // The non-empty partition
        interface.mapAndWriteData(otherMeshName, velocityName, positions, writeData);
        interface.mapAndReadData(otherMeshName, uncertaintyName, positions, dt, writeData);
      } else {
        // The empty partition
        BOOST_TEST(interface.getMeshVertexSize(otherMeshName) == 0);
        // Not possible due to the configuration
        BOOST_CHECK_THROW(interface.mapAndReadData(otherMeshName, velocityName, positions, dt, writeData), ::precice::Error);
        BOOST_CHECK_THROW(interface.mapAndWriteData(otherMeshName, uncertaintyName, positions, writeData), ::precice::Error);
        // Not possible due to the empty partition/access region
        BOOST_CHECK_THROW(interface.mapAndWriteData(otherMeshName, velocityName, positions, writeData), ::precice::Error);
        BOOST_CHECK_THROW(interface.mapAndReadData(otherMeshName, uncertaintyName, positions, dt, writeData), ::precice::Error);
      }
      interface.advance(dt);
      dt = interface.getMaxTimeStepSize();
    }
  } else {
    // Defines the mesh and reads data
    BOOST_REQUIRE(context.isNamed("SolverTwo"));
    precice::Participant interface(context.name, context.config(), context.rank, context.size);

    // Get IDs
    auto      meshName     = "MeshTwo";
    auto      velocityName = "Velocities";
    const int dim          = interface.getMeshDimensions(meshName);
    BOOST_TEST(dim == 2);
    // Define the interface
    std::vector<double> positions = context.isPrimary() ? std::vector<double>({0.0, 1.0, 0.0, 2.0}) : std::vector<double>({0.0, 3.0, 0.0, 4.0, 0.0, 5.0});

    const int        size = positions.size() / dim;
    std::vector<int> ids(size);

    interface.setMeshVertices(meshName, positions, ids);

    // Initialize the Participant
    interface.initialize();
    double dt = interface.getMaxTimeStepSize();

    // Start the time loop
    std::vector<double> readData(size);
    while (interface.isCouplingOngoing()) {

      interface.advance(dt);
      double dt = interface.getMaxTimeStepSize();
      interface.readData(meshName, velocityName, ids, dt, readData);

      // Check the received data
      const std::vector<double> expectedReadData = context.isPrimary() ? std::vector<double>({1, 2}) : expectedReadDataSecondaryRank;
    }
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Parallel
BOOST_AUTO_TEST_SUITE_END() // JustInTimeMapping

#endif // PRECICE_NO_MPI
