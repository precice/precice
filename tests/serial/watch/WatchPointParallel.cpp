#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>
#include "helpers.hpp"
#include "io/TXTTableWriter.hpp"

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
PRECICE_TEST_SETUP("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank))
BOOST_AUTO_TEST_CASE(WatchPointParallel)
{
  PRECICE_TEST();

  using Eigen::Vector2d;

  if (context.isNamed("SolverOne")) {
    precice::Participant interface(context.name, context.config(), 0, 1);

    // Set mesh
    Vector2d coord{0.0, 0.0};

    auto meshName = "MeshOne";

    int id = interface.setMeshVertex(meshName, coord);

    // Initialize, the mesh
    interface.initialize();
    double dt = interface.getMaxTimeStepSize();

    auto   dataOneID = "DataOne";
    double value     = 2.0;
    double increment = 1.0;

    while (interface.isCouplingOngoing()) {

      interface.writeData(meshName, dataOneID, {&id, 1}, {&value, 1});

      interface.advance(dt);
      double dt = interface.getMaxTimeStepSize();

      value += increment;
    }
    interface.finalize();
  } else if (context.isNamed("SolverTwo")) {

    precice::Participant interface(context.name, context.config(), 0, 1);

    // Set mesh
    Vector2d coord{0.0, 0.0};

    auto meshTwoID = "MeshTwo";

    int id = interface.setMeshVertex(meshTwoID, coord);

    // Initialize the mesh
    interface.initialize();
    double dt = interface.getMaxTimeStepSize();

    auto   dataOneID = "DataOne";
    double value;

    while (interface.isCouplingOngoing()) {

      interface.readData(meshTwoID, dataOneID, {&id, 1}, dt, {&value, 1});

      interface.advance(dt);
      double dt = interface.getMaxTimeStepSize();
    }
    interface.finalize();

    {
      std::string fileName = "precice-SolverTwo-watchpoint-WatchPoint.log";
      auto        result   = readDoublesFromTXTFile(fileName, 4);
      auto        expected = std::vector<double>{
          1.0, 0.0, 0.0, 2.0,
          2.0, 0.0, 0.0, 3.0,
          3.0, 0.0, 0.0, 4.0};
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
