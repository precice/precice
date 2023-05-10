#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(MultipleMappings)
BOOST_AUTO_TEST_CASE(MultipleReadFromMappings)
{
  PRECICE_TEST("A"_on(1_rank), "B"_on(1_rank));

  using Eigen::Vector2d;

  precice::SolverInterface interface(context.name, context.config(), context.rank, context.size);
  Vector2d                 vertex{0.0, 0.0};

  if (context.isNamed("A")) {
    auto meshNameTop    = "MeshATop";
    auto meshNameBottom = "MeshABottom";
    int  vertexIDTop    = interface.setMeshVertex(meshNameTop, vertex);
    int  vertexIDBottom = interface.setMeshVertex(meshNameBottom, vertex);
    auto dataNameTop    = "Pressure";
    auto dataNameBottom = "Pressure";

    interface.initialize();
    double dt = interface.getMaxTimeStepSize();
    interface.advance(dt);
    double pressure = -1.0;
    dt              = interface.getMaxTimeStepSize();
    interface.readData(meshNameTop, dataNameTop, {&vertexIDTop, 1}, dt, {&pressure, 1});
    BOOST_TEST(pressure == 1.0);
    pressure = -1.0;
    interface.readData(meshNameBottom, dataNameBottom, {&vertexIDBottom, 1}, dt, {&pressure, 1});
    BOOST_TEST(pressure == 1.0);
    BOOST_TEST(not interface.isCouplingOngoing());
    interface.finalize();

  } else {
    BOOST_TEST(context.isNamed("B"));
    auto meshName = "MeshB";
    int  vertexID = interface.setMeshVertex(meshName, vertex);
    auto dataName = "Pressure";

    interface.initialize();
    double pressure = 1.0;
    interface.writeData(meshName, dataName, {&vertexID, 1}, {&pressure, 1});
    double dt = interface.getMaxTimeStepSize();
    interface.advance(dt);
    BOOST_TEST(not interface.isCouplingOngoing());
    interface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // MultipleMappings

#endif // PRECICE_NO_MPI
