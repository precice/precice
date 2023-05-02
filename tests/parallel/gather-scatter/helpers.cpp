#ifndef PRECICE_NO_MPI

#include "helpers.hpp"

#include "precice/SolverInterface.hpp"
#include "testing/Testing.hpp"

// In order to test enforced gather scatter communication with an empty primary rank (see below)
void runTestEnforceGatherScatter(std::vector<double> primaryPartition, const TestContext &context)
{
  if (context.isNamed("ParallelSolver")) {
    // Get mesh and data IDs
    precice::SolverInterface interface(context.name, context.config(), context.rank, context.size);
    auto                     meshName      = "ParallelMesh";
    auto                     writeDataName = "MyData1";
    auto                     readDataName  = "MyData2";
    const int                dim           = interface.getMeshDimensions(meshName);
    BOOST_TEST(dim == 2);

    // Set coordinates, primary according to input argument
    const std::vector<double> coordinates = context.isPrimary() ? primaryPartition : std::vector<double>{0.0, 0.5, 0.0, 3.5, 0.0, 5.0};
    const unsigned int        size        = coordinates.size() / dim;
    std::vector<int>          ids(size, 0);

    // Set mesh vertices
    interface.setMeshVertices(meshName, size, coordinates.data(), ids.data());

    // Initialize the solverinterface
    interface.initialize();
    double dt = interface.getMaxTimeStepSize();

    // Create some dummy writeData
    std::vector<double> writeData;
    for (unsigned int i = 0; i < size; ++i) {
      writeData.emplace_back(i + 1);
    }

    // Allocate memory for readData
    std::vector<double> readData(size);
    while (interface.isCouplingOngoing()) {
      // Write data, advance the solverinterface and readData
      interface.writeBlockScalarData(meshName, writeDataName, size,
                                     ids.data(), writeData.data());

      interface.advance(dt);
      double dt = interface.getMaxTimeStepSize();
      interface.readBlockScalarData(meshName, readDataName, size,
                                    ids.data(), dt, readData.data());
      // The received data on the secondary rank is always the same
      if (!context.isPrimary()) {
        BOOST_TEST(readData == std::vector<double>({3.4, 5.7, 4.0}));
      }
    }
  } else {
    // The serial participant
    BOOST_REQUIRE(context.isNamed("SerialSolver"));
    precice::SolverInterface interface(context.name, context.config(), context.rank, context.size);
    // Get IDs
    auto      meshName      = "SerialMesh";
    auto      writeDataName = "MyData2";
    auto      readDataName  = "MyData1";
    const int dim           = interface.getMeshDimensions(meshName);
    BOOST_TEST(interface.getMeshDimensions(meshName) == 2);

    // Define the interface
    const std::vector<double> coordinates{0.0, 0.5, 0.0, 3.5, 0.0, 5.0};
    const unsigned int        size = coordinates.size() / dim;
    std::vector<int>          ids(size);

    // Set vertices
    interface.setMeshVertices(meshName, size, coordinates.data(), ids.data());

    // Initialize the solverinterface
    interface.initialize();
    double dt = interface.getMaxTimeStepSize();

    // Somce arbitrary write data
    std::vector<double> writeData{3.4, 5.7, 4.0};
    std::vector<double> readData(size);

    // Start the time loop
    while (interface.isCouplingOngoing()) {
      // Write data, advance solverinterface and read data
      interface.writeBlockScalarData(meshName, writeDataName, size,
                                     ids.data(), writeData.data());
      interface.advance(dt);
      double dt = interface.getMaxTimeStepSize();
      interface.readBlockScalarData(meshName, readDataName, size,
                                    ids.data(), dt, readData.data());
      // The received data is always the same
      if (!context.isPrimary()) {
        BOOST_TEST(readData == std::vector<double>({1, 2, 3}));
      }
    }
  }
}

#endif
