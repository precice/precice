#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
#include <vector>
#include "helpers.hpp"
#include "io/TXTTableWriter.hpp"

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_CASE(WatchIntegralScaleAndNoScale)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  using Eigen::Vector2d;

  if (context.isNamed("SolverOne")) {
    precice::SolverInterface interface(context.name, context.config(), 0, 1);

    // Set mesh
    Vector2d coordA{0.0, 0.0};
    Vector2d coordB{1.0, 0.0};
    Vector2d coordC{1.0, 2.0};

    auto meshName = "MeshOne";

    int idA = interface.setMeshVertex(meshName, coordA.data());
    int idB = interface.setMeshVertex(meshName, coordB.data());
    int idC = interface.setMeshVertex(meshName, coordC.data());

    interface.setMeshEdge(meshName, idA, idB);
    interface.setMeshEdge(meshName, idB, idC);

    // Initialize, the mesh
    double dt = interface.initialize();

    auto   dataOneID = "DataOne";
    double valueA    = 1.0;
    double valueB    = 2.0;
    double valueC    = 3.0;

    double increment = 1.0;

    while (interface.isCouplingOngoing()) {

      interface.writeScalarData(meshName, dataOneID, idA, valueA);
      interface.writeScalarData(meshName, dataOneID, idB, valueB);
      interface.writeScalarData(meshName, dataOneID, idC, valueC);

      dt = interface.advance(dt);

      valueA += increment;
      valueB += increment;
      valueC += increment;
    }
    interface.finalize();
  } else if (context.isNamed("SolverTwo")) {

    precice::SolverInterface interface(context.name, context.config(), 0, 1);

    // Set mesh
    Vector2d coordA{0.0, 0.0};
    Vector2d coordB{1.0, 0.0};
    Vector2d coordC{1.0, 2.0};

    auto meshTwoID = "MeshTwo";

    int idA = interface.setMeshVertex(meshTwoID, coordA.data());
    int idB = interface.setMeshVertex(meshTwoID, coordB.data());
    int idC = interface.setMeshVertex(meshTwoID, coordC.data());

    interface.setMeshEdge(meshTwoID, idA, idB);
    interface.setMeshEdge(meshTwoID, idB, idC);

    // Initialize the mesh
    double dt = interface.initialize();

    auto   dataOneID = "DataOne";
    double valueA, valueB, valueC;

    while (interface.isCouplingOngoing()) {

      interface.readScalarData(meshTwoID, dataOneID, idA, dt, valueA);
      interface.readScalarData(meshTwoID, dataOneID, idB, dt, valueB);
      interface.readScalarData(meshTwoID, dataOneID, idC, dt, valueC);

      dt = interface.advance(dt);
    }
    interface.finalize();

    {
      std::string fileName = "precice-SolverTwo-watchintegral-WatchIntegral.log";
      auto        result   = readDoublesFromTXTFile(fileName, 3);
      auto        expected = std::vector<double>{
          1.0, 9.5, 3.0,
          2.0, 12.5, 3.0,
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
          1.0, 9.0, 3.0,
          2.0, 12.0, 3.0,
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
