#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(MultipleMappings)
PRECICE_TEST_SETUP("A"_on(1_rank), "B"_on(1_rank))
BOOST_AUTO_TEST_CASE(MultipleWriteFromMappingsAndData)
{
  PRECICE_TEST();

  using Eigen::Vector2d;

  precice::Participant interface(context.name, context.config(), context.rank, context.size);
  Vector2d             vertex1{0.0, 0.0};
  Vector2d             vertex2{2.0, 0.0};
  Vector2d             vertex3{4.0, 0.0};

  if (context.isNamed("A")) {
    auto meshNameTop     = "MeshATop";
    auto meshNameBottom  = "MeshABottom";
    int  vertexIDTop     = interface.setMeshVertex(meshNameTop, vertex1);
    int  vertexIDBottom  = interface.setMeshVertex(meshNameBottom, vertex3);
    auto dataNameTopP    = "Pressure";
    auto dataNameBottomP = "Pressure";
    auto dataNameTopT    = "Temperature";
    auto dataNameBottomT = "Temperature";

    interface.initialize();
    double dt = interface.getMaxTimeStepSize();
    interface.advance(dt);
    double pressure    = -1.0;
    double temperature = -1.0;
    dt                 = interface.getMaxTimeStepSize();
    interface.readData(meshNameTop, dataNameTopP, {&vertexIDTop, 1}, dt, {&pressure, 1});
    interface.readData(meshNameTop, dataNameTopT, {&vertexIDTop, 1}, dt, {&temperature, 1});
    BOOST_TEST(pressure == 1.0);
    BOOST_TEST(temperature == 331);
    pressure    = -1.0;
    temperature = -1.0;
    interface.readData(meshNameBottom, dataNameBottomP, {&vertexIDBottom, 1}, dt, {&pressure, 1});
    interface.readData(meshNameBottom, dataNameBottomT, {&vertexIDBottom, 1}, dt, {&temperature, 1});
    BOOST_TEST(temperature == 273.15);
    BOOST_TEST(pressure == 5.0);
    BOOST_TEST(not interface.isCouplingOngoing());
    interface.finalize();

  } else {
    BOOST_TEST(context.isNamed("B"));
    auto meshName  = "MeshB";
    int  vertexID1 = interface.setMeshVertex(meshName, vertex1);
    int  vertexID2 = interface.setMeshVertex(meshName, vertex2);
    int  vertexID3 = interface.setMeshVertex(meshName, vertex3);
    auto dataNameP = "Pressure";
    auto dataNameT = "Temperature";

    interface.initialize();
    double pressure    = 1.0;
    double temperature = 331;
    interface.writeData(meshName, dataNameP, {&vertexID1, 1}, {&pressure, 1});
    interface.writeData(meshName, dataNameT, {&vertexID1, 1}, {&temperature, 1});
    pressure    = 4.0;
    temperature = 335;
    interface.writeData(meshName, dataNameP, {&vertexID2, 1}, {&pressure, 1});
    interface.writeData(meshName, dataNameT, {&vertexID2, 1}, {&temperature, 1});
    pressure    = 5.0;
    temperature = 273.15;
    interface.writeData(meshName, dataNameP, {&vertexID3, 1}, {&pressure, 1});
    interface.writeData(meshName, dataNameT, {&vertexID3, 1}, {&temperature, 1});
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
