#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(PreciceTests)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(MultipleMappings)
BOOST_AUTO_TEST_CASE(MultipleToMappings)
{
  PRECICE_TEST("A"_on(1_rank), "B"_on(1_rank));

  using Eigen::Vector2d;
  using namespace precice::constants;

  precice::SolverInterface interface(context.name, context.config(), context.rank, context.size);
  Vector2d                 vertex{0.0, 0.0};

  if (context.isNamed("A")) {
    const precice::MeshID meshIDTop      = interface.getMeshID("MeshATop");
    const precice::MeshID meshIDBottom   = interface.getMeshID("MeshABottom");
    int                   vertexIDTop    = interface.setMeshVertex(meshIDTop, vertex.data());
    int                   vertexIDBottom = interface.setMeshVertex(meshIDBottom, vertex.data());
    int                   dataIDTop      = interface.getDataID("DisplacementTop", meshIDTop);
    int                   dataIDBottom   = interface.getDataID("DisplacementBottom", meshIDBottom);

    double dt              = interface.initialize();
    double displacementTop = 1.0;
    interface.writeScalarData(dataIDTop, vertexIDTop, displacementTop);
    double displacementBottom = 2.0;
    interface.writeScalarData(dataIDBottom, vertexIDBottom, displacementBottom);
    interface.advance(dt);
    BOOST_TEST(not interface.isCouplingOngoing());
    interface.finalize();

  } else {
    BOOST_TEST(context.isNamed("B"));
    const precice::MeshID meshID   = interface.getMeshID("MeshB");
    int                   vertexID = interface.setMeshVertex(meshID, vertex.data());
    int                   dataID   = interface.getDataID("DisplacementSum", meshID);

    double dt = interface.initialize();
    interface.advance(dt);
    double displacement = -1.0;
    interface.readScalarData(dataID, vertexID, displacement);
    BOOST_TEST(displacement == 3.0);
    BOOST_TEST(not interface.isCouplingOngoing());
    interface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END() // PreciceTests
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // MultipleMappings

#endif // PRECICE_NO_MPI
