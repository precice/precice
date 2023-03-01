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
    auto meshIDTop      = "MeshATop";
    auto meshIDBottom   = "MeshABottom";
    int  vertexIDTop    = interface.setMeshVertex(meshIDTop, vertex1.data());
    int  vertexIDBottom = interface.setMeshVertex(meshIDBottom, vertex3.data());
    auto dataIDTopP     = "Pressure";    //  meshIDTop
    auto dataIDBottomP  = "Pressure";    //  meshIDBottom
    auto dataIDTopT     = "Temperature"; //  meshIDTop
    auto dataIDBottomT  = "Temperature"; //  meshIDBottom

    double dt = interface.initialize();
    interface.advance(dt);
    double pressure    = -1.0;
    double temperature = -1.0;
    interface.readScalarData(meshIDTop, dataIDTopP, vertexIDTop, pressure);
    interface.readScalarData(meshIDTop, dataIDTopT, vertexIDTop, temperature);
    BOOST_TEST(pressure == 1.0);
    BOOST_TEST(temperature == 331);
    pressure    = -1.0;
    temperature = -1.0;
    interface.readScalarData(meshIDBottom, dataIDBottomP, vertexIDBottom, pressure);
    interface.readScalarData(meshIDBottom, dataIDBottomT, vertexIDBottom, temperature);
    BOOST_TEST(temperature == 273.15);
    BOOST_TEST(pressure == 5.0);
    BOOST_TEST(not interface.isCouplingOngoing());
    interface.finalize();

  } else {
    BOOST_TEST(context.isNamed("B"));
    auto meshID    = "MeshB";
    int  vertexID1 = interface.setMeshVertex(meshID, vertex1.data());
    int  vertexID2 = interface.setMeshVertex(meshID, vertex2.data());
    int  vertexID3 = interface.setMeshVertex(meshID, vertex3.data());
    auto dataIDP   = "Pressure";    //  meshID
    auto dataIDT   = "Temperature"; //  meshID

    double dt          = interface.initialize();
    double pressure    = 1.0;
    double temperature = 331;
    interface.writeScalarData(meshID, dataIDP, vertexID1, pressure);
    interface.writeScalarData(meshID, dataIDT, vertexID1, temperature);
    pressure    = 4.0;
    temperature = 335;
    interface.writeScalarData(meshID, dataIDP, vertexID2, pressure);
    interface.writeScalarData(meshID, dataIDT, vertexID2, temperature);
    pressure    = 5.0;
    temperature = 273.15;
    interface.writeScalarData(meshID, dataIDP, vertexID3, pressure);
    interface.writeScalarData(meshID, dataIDT, vertexID3, temperature);
    interface.advance(dt);
    BOOST_TEST(not interface.isCouplingOngoing());
    interface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // MultipleMappings

#endif // PRECICE_NO_MPI
