#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>
#include "helpers.hpp"
#include "io/TXTTableWriter.hpp"

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
PRECICE_TEST_SETUP("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank))
BOOST_AUTO_TEST_CASE(WatchIntegralScaleAndNoScaleParallel)
{
  PRECICE_TEST();

  if (context.isNamed("SolverOne")) {
    precice::Participant interface(context.name, context.config(), 0, 1);

    // Set mesh
    std::vector<double> coords{
        0.0, 0.0,
        1.0, 0.0,
        1.0, 2.0};

    auto meshName = "MeshOne";

    const int        nVertices = 3;
    std::vector<int> ids(nVertices);
    interface.setMeshVertices(meshName, coords, ids);

    interface.setMeshEdge(meshName, ids[0], ids[1]);
    interface.setMeshEdge(meshName, ids[1], ids[2]);

    // Initialize, the mesh
    interface.initialize();
    double dt = interface.getMaxTimeStepSize();

    auto                dataOneID = "DataOne";
    std::vector<double> values{1.0, 2.0, 3.0};

    double increment = 1.0;

    while (interface.isCouplingOngoing()) {

      interface.writeData(meshName, dataOneID, ids, values);

      interface.advance(dt);

      for (int i = 0; i < nVertices; i++) {
        values[i] += increment;
      }
    }
    interface.finalize();
  } else if (context.isNamed("SolverTwo")) {

    precice::Participant interface(context.name, context.config(), 0, 1);

    // Set mesh
    std::vector<double> coords{
        0.0, 0.0,
        1.0, 0.0,
        1.0, 2.0};
    auto meshName = "MeshTwo";

    const int        nVertices = 3;
    std::vector<int> ids(nVertices);
    interface.setMeshVertices(meshName, coords, ids);

    interface.setMeshEdge(meshName, ids[0], ids[1]);
    interface.setMeshEdge(meshName, ids[1], ids[2]);

    // Initialize the mesh
    interface.initialize();
    while (interface.isCouplingOngoing()) {
      interface.advance(interface.getMaxTimeStepSize());
    }
    interface.finalize();

    {
      std::string fileName = "precice-SolverTwo-watchintegral-WatchIntegral.log";
      auto        result   = readDoublesFromTXTFile(fileName, 3);
      auto        expected = std::vector<double>{
          // Time  DataOne  SurfaceArea
          0.0, 0.0, 3.0,
          1.0, 6.5, 3.0,
          2.0, 9.5, 3.0,
          3.0, 12.5, 3.0};
      BOOST_TEST(result.size() == expected.size());
      for (size_t i = 0; i < result.size(); ++i) {
        BOOST_TEST_CONTEXT("entry index: " << i)
        {
          using testing::equals;
          BOOST_TEST(equals(result.at(i), expected.at(i)));
        }
      }
    }

    {
      std::string fileName = "precice-SolverTwo-watchintegral-WatchIntegralNoScale.log";
      auto        result   = readDoublesFromTXTFile(fileName, 3);
      auto        expected = std::vector<double>{
          // Time  DataOne  SurfaceArea
          0.0, 0.0, 3.0,
          1.0, 6.0, 3.0,
          2.0, 9.0, 3.0,
          3.0, 12.0, 3.0};
      BOOST_TEST(result.size() == expected.size());
      for (size_t i = 0; i < result.size(); ++i) {
        BOOST_TEST_CONTEXT("entry index: " << i)
        {
          using testing::equals;
          BOOST_TEST(equals(result.at(i), expected.at(i)));
        }
      }
    }
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial

#endif // PRECICE_NO_MPI
