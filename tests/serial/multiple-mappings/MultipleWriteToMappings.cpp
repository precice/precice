#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(MultipleMappings)
PRECICE_TEST_SETUP("A"_on(1_rank), "B"_on(1_rank))
BOOST_AUTO_TEST_CASE(MultipleWriteToMappings)
{
  PRECICE_TEST();

  using Eigen::Vector2d;

  precice::Participant interface(context.name, context.config(), context.rank, context.size);
  Vector2d             vertex{0.0, 0.0};

  if (context.isNamed("A")) {
    auto meshNameTop    = "MeshATop";
    auto meshNameBottom = "MeshABottom";
    int  vertexIDTop    = interface.setMeshVertex(meshNameTop, vertex);
    int  vertexIDBottom = interface.setMeshVertex(meshNameBottom, vertex);
    auto dataNameTop    = "DisplacementTop";
    auto dataNameBottom = "DisplacementBottom";

    interface.initialize();
    double displacementTop = 1.0;
    interface.writeData(meshNameTop, dataNameTop, {&vertexIDTop, 1}, {&displacementTop, 1});
    double displacementBottom = 2.0;
    interface.writeData(meshNameBottom, dataNameBottom, {&vertexIDBottom, 1}, {&displacementBottom, 1});
    double dt = interface.getMaxTimeStepSize();
    interface.advance(dt);
    BOOST_TEST(not interface.isCouplingOngoing());
    interface.finalize();

  } else {
    BOOST_TEST(context.isNamed("B"));
    auto meshName = "MeshB";
    int  vertexID = interface.setMeshVertex(meshName, vertex);
    auto dataName = "DisplacementSum";

    interface.initialize();
    double dt = interface.getMaxTimeStepSize();
    interface.advance(dt);
    double displacement = -1.0;
    dt                  = interface.getMaxTimeStepSize();
    interface.readData(meshName, dataName, {&vertexID, 1}, dt, {&displacement, 1});
    BOOST_TEST(displacement == 3.0);
    BOOST_TEST(not interface.isCouplingOngoing());
    interface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // MultipleMappings

#endif // PRECICE_NO_MPI
