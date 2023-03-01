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
    auto meshIDTop      = "MeshATop";
    auto meshIDBottom   = "MeshABottom";
    int  vertexIDTop    = interface.setMeshVertex(meshIDTop, vertex1.data());
    int  vertexIDBottom = interface.setMeshVertex(meshIDBottom, vertex3.data());
    auto dataIDTop      = "Pressure"; //  meshIDTop
    auto dataIDBottom   = "Pressure"; //  meshIDBottom

    double dt = interface.initialize();
    interface.advance(dt);
    double pressure = -1.0;
    interface.readScalarData(dataIDTop, vertexIDTop, pressure);
    BOOST_TEST(pressure == 1.0);
    pressure = -1.0;
    interface.readScalarData(dataIDBottom, vertexIDBottom, pressure);
    BOOST_TEST(pressure == 5.0);
    BOOST_TEST(not interface.isCouplingOngoing());
    interface.finalize();

  } else {
    BOOST_TEST(context.isNamed("B"));
    auto meshID    = "MeshB";
    int  vertexID1 = interface.setMeshVertex(meshID, vertex1.data());
    int  vertexID2 = interface.setMeshVertex(meshID, vertex2.data());
    int  vertexID3 = interface.setMeshVertex(meshID, vertex3.data());
    auto dataID    = "Pressure"; //  meshID

    double dt       = interface.initialize();
    double pressure = 1.0;
    interface.writeScalarData(dataID, vertexID1, pressure);
    pressure = 4.0;
    interface.writeScalarData(dataID, vertexID2, pressure);
    pressure = 5.0;
    interface.writeScalarData(dataID, vertexID3, pressure);
    interface.advance(dt);
    BOOST_TEST(not interface.isCouplingOngoing());
    interface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // MultipleMappings

#endif // PRECICE_NO_MPI
