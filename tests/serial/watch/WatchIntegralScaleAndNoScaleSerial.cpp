#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>
#include "helpers.hpp"
#include "io/TXTTableWriter.hpp"

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
PRECICE_TEST_SETUP("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank))
BOOST_AUTO_TEST_CASE(WatchIntegralScaleAndNoScaleSerial)
{
  PRECICE_TEST();

  using Eigen::Vector2d;

  if (context.isNamed("SolverOne")) {
    precice::Participant interface(context.name, context.config(), 0, 1);

    // Set mesh
    Vector2d coordA{0.0, 0.0};
    Vector2d coordB{1.0, 0.0};
    Vector2d coordC{1.0, 2.0};

    auto meshName = "MeshOne";

    int idA = interface.setMeshVertex(meshName, coordA);
    int idB = interface.setMeshVertex(meshName, coordB);
    int idC = interface.setMeshVertex(meshName, coordC);

    interface.setMeshEdge(meshName, idA, idB);
    interface.setMeshEdge(meshName, idB, idC);

    // Initialize, the mesh
    interface.initialize();
    double dt = interface.getMaxTimeStepSize();

    auto   dataOneID = "DataOne";
    double valueA    = 1.0;
    double valueB    = 2.0;
    double valueC    = 3.0;

    double increment = 1.0;

    while (interface.isCouplingOngoing()) {

      interface.writeData(meshName, dataOneID, {&idA, 1}, {&valueA, 1});
      interface.writeData(meshName, dataOneID, {&idB, 1}, {&valueB, 1});
      interface.writeData(meshName, dataOneID, {&idC, 1}, {&valueC, 1});

      interface.advance(dt);
      double dt = interface.getMaxTimeStepSize();

      valueA += increment;
      valueB += increment;
      valueC += increment;
    }
    interface.finalize();
  } else if (context.isNamed("SolverTwo")) {

    precice::Participant interface(context.name, context.config(), 0, 1);

    // Set mesh
    Vector2d coordA{0.0, 0.0};
    Vector2d coordB{1.0, 0.0};
    Vector2d coordC{1.0, 2.0};

    auto meshTwoID = "MeshTwo";

    int idA = interface.setMeshVertex(meshTwoID, coordA);
    int idB = interface.setMeshVertex(meshTwoID, coordB);
    int idC = interface.setMeshVertex(meshTwoID, coordC);

    interface.setMeshEdge(meshTwoID, idA, idB);
    interface.setMeshEdge(meshTwoID, idB, idC);

    // Initialize the mesh
    interface.initialize();
    double dt = interface.getMaxTimeStepSize();

    auto   dataOneID = "DataOne";
    double valueA, valueB, valueC;

    while (interface.isCouplingOngoing()) {

      interface.readData(meshTwoID, dataOneID, {&idA, 1}, dt, {&valueA, 1});
      interface.readData(meshTwoID, dataOneID, {&idB, 1}, dt, {&valueB, 1});
      interface.readData(meshTwoID, dataOneID, {&idC, 1}, dt, {&valueC, 1});

      interface.advance(dt);
      double dt = interface.getMaxTimeStepSize();
    }
    interface.finalize();

    {
      std::string fileName = "precice-SolverTwo-watchintegral-WatchIntegral.log";
      auto        result   = readDoublesFromTXTFile(fileName, 3);
      auto        expected = std::vector<double>{
          // Time  DataOne  SurfaceArea
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
