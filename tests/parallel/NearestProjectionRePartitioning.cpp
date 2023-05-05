#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <boost/test/data/test_case.hpp>
#include <precice/SolverInterface.hpp>
#include <precice/impl/SolverInterfaceImpl.hpp>
#include <vector>

/// This testcase is based on a bug documented in issue #371
BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Parallel)
BOOST_DATA_TEST_CASE(NearestProjectionRePartitioning,
                     boost::unit_test::data::make({true, false}),
                     useBulkFunctions)
{
  PRECICE_TEST("FluidSolver"_on(3_ranks), "SolidSolver"_on(1_rank));

  if (context.isNamed("FluidSolver")) {
    precice::SolverInterface interface(context.name, context.config(), context.rank, context.size);

    if (context.isPrimary()) {
      interface.initialize();
      interface.advance(1.0);
      interface.finalize();
    } else {
      auto      meshName   = "CellCenters";
      const int dimensions = interface.getMeshDimensions(meshName);
      BOOST_REQUIRE(dimensions == 3);

      const int                 numberOfVertices = 65;
      const double              yCoord           = 0.0;
      const double              zCoord           = 0.005;
      const std::vector<double> positions{
          0.00124795, yCoord, zCoord,
          0.00375646, yCoord, zCoord,
          0.00629033, yCoord, zCoord,
          0.00884982, yCoord, zCoord,
          0.0114352, yCoord, zCoord,
          0.0140467, yCoord, zCoord,
          0.0166846, yCoord, zCoord,
          0.0193492, yCoord, zCoord,
          0.0220407, yCoord, zCoord,
          0.0247594, yCoord, zCoord,
          0.0275056, yCoord, zCoord,
          0.0302796, yCoord, zCoord,
          0.0330816, yCoord, zCoord,
          0.0359119, yCoord, zCoord,
          0.0387709, yCoord, zCoord,
          0.0416588, yCoord, zCoord,
          0.0445758, yCoord, zCoord,
          0.0475224, yCoord, zCoord,
          0.0504987, yCoord, zCoord,
          0.0535051, yCoord, zCoord,
          0.0565419, yCoord, zCoord,
          0.0596095, yCoord, zCoord,
          0.062708, yCoord, zCoord,
          0.0658378, yCoord, zCoord,
          0.0689993, yCoord, zCoord,
          0.0721928, yCoord, zCoord,
          0.0754186, yCoord, zCoord,
          0.0786769, yCoord, zCoord,
          0.0819682, yCoord, zCoord,
          0.0852928, yCoord, zCoord,
          0.088651, yCoord, zCoord,
          0.0920431, yCoord, zCoord,
          0.0954695, yCoord, zCoord,
          0.0989306, yCoord, zCoord,
          0.102427, yCoord, zCoord,
          0.105958, yCoord, zCoord,
          0.109525, yCoord, zCoord,
          0.113128, yCoord, zCoord,
          0.116768, yCoord, zCoord,
          0.120444, yCoord, zCoord,
          0.124158, yCoord, zCoord,
          0.127909, yCoord, zCoord,
          0.131698, yCoord, zCoord,
          0.135525, yCoord, zCoord,
          0.139391, yCoord, zCoord,
          0.143296, yCoord, zCoord,
          0.147241, yCoord, zCoord,
          0.151226, yCoord, zCoord,
          0.15525, yCoord, zCoord,
          0.159316, yCoord, zCoord,
          0.163422, yCoord, zCoord,
          0.16757, yCoord, zCoord,
          0.17176, yCoord, zCoord,
          0.175993, yCoord, zCoord,
          0.180268, yCoord, zCoord,
          0.184586, yCoord, zCoord,
          0.188948, yCoord, zCoord,
          0.193354, yCoord, zCoord,
          0.197805, yCoord, zCoord,
          0.202301, yCoord, zCoord,
          0.206842, yCoord, zCoord,
          0.211429, yCoord, zCoord,
          0.216062, yCoord, zCoord,
          0.220742, yCoord, zCoord,
          0.22547, yCoord, zCoord};
      BOOST_TEST(numberOfVertices * dimensions == positions.size());
      std::vector<int> vertexIDs(numberOfVertices);
      interface.setMeshVertices(meshName, positions, vertexIDs);
      interface.initialize();
      BOOST_TEST(precice::testing::WhiteboxAccessor::impl(interface).mesh("Nodes").triangles().size() == 15);
      interface.advance(1.0);
      interface.finalize();
    }
  } else {
    BOOST_TEST(context.isNamed("SolidSolver"));
    precice::SolverInterface interface(context.name, context.config(), context.rank, context.size);
    auto                     meshName   = "Nodes";
    const int                dimensions = interface.getMeshDimensions(meshName);
    BOOST_REQUIRE(dimensions == 3);
    const int                 numberOfVertices = 34;
    const double              yCoord           = 0.0;
    const double              zCoord1          = 0.0;
    const double              zCoord2          = 0.01;
    const std::vector<double> positions{
        0.0, yCoord, zCoord2,
        0.0, yCoord, zCoord1,
        0.03125, yCoord, zCoord2,
        0.03125, yCoord, zCoord1,
        0.0625, yCoord, zCoord2,
        0.0625, yCoord, zCoord1,
        0.09375, yCoord, zCoord2,
        0.09375, yCoord, zCoord1,
        0.125, yCoord, zCoord2,
        0.125, yCoord, zCoord1,
        0.15625, yCoord, zCoord2,
        0.15625, yCoord, zCoord1,
        0.1875, yCoord, zCoord2,
        0.1875, yCoord, zCoord1,
        0.21875, yCoord, zCoord2,
        0.21875, yCoord, zCoord1,
        0.25, yCoord, zCoord2,
        0.25, yCoord, zCoord1,
        0.28125, yCoord, zCoord2,
        0.28125, yCoord, zCoord1,
        0.3125, yCoord, zCoord2,
        0.3125, yCoord, zCoord1,
        0.34375, yCoord, zCoord2,
        0.34375, yCoord, zCoord1,
        0.375, yCoord, zCoord2,
        0.375, yCoord, zCoord1,
        0.40625, yCoord, zCoord2,
        0.40625, yCoord, zCoord1,
        0.4375, yCoord, zCoord2,
        0.4375, yCoord, zCoord1,
        0.46875, yCoord, zCoord2,
        0.46875, yCoord, zCoord1,
        0.5, yCoord, zCoord2,
        0.5, yCoord, zCoord1};
    BOOST_TEST(numberOfVertices * dimensions == positions.size());
    std::vector<int> vertexIDs(numberOfVertices);
    interface.setMeshVertices(meshName, positions, vertexIDs);

    const int numberOfCells = numberOfVertices / 2 - 1;

    if (useBulkFunctions) {
      std::vector<int> ids;
      for (int i = 0; i < numberOfCells; i++) {
        // left-diag-bottom
        ids.push_back(vertexIDs.at(i * 2));
        ids.push_back(vertexIDs.at(i * 2 + 1));
        ids.push_back(vertexIDs.at(i * 2 + 3));

        // top-diag-right
        ids.push_back(vertexIDs.at(i * 2));
        ids.push_back(vertexIDs.at(i * 2 + 2));
        ids.push_back(vertexIDs.at(i * 2 + 3));
      }
      interface.setMeshTriangles(meshName, ids);
    } else {
      for (int i = 0; i < numberOfCells; i++) {
        // left-diag-bottom
        interface.setMeshTriangle(meshName, vertexIDs.at(i * 2), vertexIDs.at(i * 2 + 1), vertexIDs.at(i * 2 + 3));

        // top-diag-right
        interface.setMeshTriangle(meshName, vertexIDs.at(i * 2), vertexIDs.at(i * 2 + 2), vertexIDs.at(i * 2 + 3));
      }
    }

    interface.initialize();
    interface.advance(1.0);
    interface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Parallel

#endif // PRECICE_NO_MPI
