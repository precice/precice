#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(MultipleMappings)
BOOST_AUTO_TEST_CASE(MultipleWriteFromMappingsAndData)
{
  PRECICE_TEST("A"_on(1_rank), "B"_on(1_rank));

  using Eigen::Vector2d;

  precice::SolverInterface interface(context.name, context.config(), context.rank, context.size);
  Vector2d                 vertex1{0.0, 0.0};
  Vector2d                 vertex2{2.0, 0.0};
  Vector2d                 vertex3{4.0, 0.0};

  if (context.isNamed("A")) {
    auto meshNameTop     = "MeshATop";
    auto meshNameBottom  = "MeshABottom";
    int  vertexIDTop     = interface.setMeshVertex(meshNameTop, vertex1.data());
    int  vertexIDBottom  = interface.setMeshVertex(meshNameBottom, vertex3.data());
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
    interface.readScalarData(meshNameTop, dataNameTopP, vertexIDTop, dt, pressure);
    interface.readScalarData(meshNameTop, dataNameTopT, vertexIDTop, dt, temperature);
    BOOST_TEST(pressure == 1.0);
    BOOST_TEST(temperature == 331);
    pressure    = -1.0;
    temperature = -1.0;
    interface.readScalarData(meshNameBottom, dataNameBottomP, vertexIDBottom, dt, pressure);
    interface.readScalarData(meshNameBottom, dataNameBottomT, vertexIDBottom, dt, temperature);
    BOOST_TEST(temperature == 273.15);
    BOOST_TEST(pressure == 5.0);
    BOOST_TEST(not interface.isCouplingOngoing());
    interface.finalize();

  } else {
    BOOST_TEST(context.isNamed("B"));
    auto meshName  = "MeshB";
    int  vertexID1 = interface.setMeshVertex(meshName, vertex1.data());
    int  vertexID2 = interface.setMeshVertex(meshName, vertex2.data());
    int  vertexID3 = interface.setMeshVertex(meshName, vertex3.data());
    auto dataNameP = "Pressure";
    auto dataNameT = "Temperature";

    interface.initialize();
    double pressure    = 1.0;
    double temperature = 331;
    interface.writeScalarData(meshName, dataNameP, vertexID1, pressure);
    interface.writeScalarData(meshName, dataNameT, vertexID1, temperature);
    pressure    = 4.0;
    temperature = 335;
    interface.writeScalarData(meshName, dataNameP, vertexID2, pressure);
    interface.writeScalarData(meshName, dataNameT, vertexID2, temperature);
    pressure    = 5.0;
    temperature = 273.15;
    interface.writeScalarData(meshName, dataNameP, vertexID3, pressure);
    interface.writeScalarData(meshName, dataNameT, vertexID3, temperature);
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
