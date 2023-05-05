#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(MultipleMappings)
BOOST_AUTO_TEST_CASE(MultipleWriteFromMappings)
{
  PRECICE_TEST("A"_on(1_rank), "B"_on(1_rank));

  using Eigen::Vector2d;

  precice::SolverInterface interface(context.name, context.config(), context.rank, context.size);
  Vector2d                 vertex1{0.0, 0.0};
  Vector2d                 vertex2{2.0, 0.0};
  Vector2d                 vertex3{4.0, 0.0};

  if (context.isNamed("A")) {
    auto meshNameTop    = "MeshATop";
    auto meshNameBottom = "MeshABottom";
    int  vertexIDTop    = interface.setMeshVertex(meshNameTop, vertex1.data());
    int  vertexIDBottom = interface.setMeshVertex(meshNameBottom, vertex3.data());
    auto dataNameTop    = "Pressure";
    auto dataNameBottom = "Pressure";

    interface.initialize();
    double dt = interface.getMaxTimeStepSize();
    interface.advance(dt);
    dt              = interface.getMaxTimeStepSize();
    double pressure = -1.0;
    interface.readData(meshNameTop, dataNameTop, {&vertexIDTop, 1}, dt, {&pressure, 1});
    BOOST_TEST(pressure == 1.0);
    pressure = -1.0;
    interface.readData(meshNameBottom, dataNameBottom, {&vertexIDBottom, 1}, dt, {&pressure, 1});
    BOOST_TEST(pressure == 5.0);
    BOOST_TEST(not interface.isCouplingOngoing());
    interface.finalize();

  } else {
    BOOST_TEST(context.isNamed("B"));
    auto meshName  = "MeshB";
    int  vertexID1 = interface.setMeshVertex(meshName, vertex1.data());
    int  vertexID2 = interface.setMeshVertex(meshName, vertex2.data());
    int  vertexID3 = interface.setMeshVertex(meshName, vertex3.data());
    auto dataName  = "Pressure";

    interface.initialize();
    double pressure = 1.0;
    interface.writeData(meshName, dataName, {&vertexID1, 1}, {&pressure, 1});
    pressure = 4.0;
    interface.writeData(meshName, dataName, {&vertexID2, 1}, {&pressure, 1});
    pressure = 5.0;
    interface.writeData(meshName, dataName, {&vertexID3, 1}, {&pressure, 1});
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
