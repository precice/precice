#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <boost/test/tools/detail/per_element_manip.hpp>
#include <precice/SolverInterface.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(Compositional)
BOOST_AUTO_TEST_CASE(OneActivatedMuscle)
{
  PRECICE_TEST("M1SM"_on(1_rank), "M2SM"_on(1_rank), "M1"_on(1_rank), "Tendon"_on(1_rank));

  std::cout << "Before constructor" << std::endl;
  precice::SolverInterface interface(context.name, context.config(), context.rank, context.size);

  const std::vector<double> surfaceCoords{1, 0, 2, 0};
  const std::vector<double> signalCoords{0, 0};

  int              activationDataID, stretchDataID, displacement1DataID, traction1DataID, displacement2DataID, traction2DataID, crossStretchDataID;
  std::vector<int> surface1VertexIDs(2);
  std::vector<int> surface2VertexIDs(2);
  std::vector<int> activationVertexIDs(1);
  std::vector<int> stretchVertexIDs(1);
  std::vector<int> crossStretchVertexIDs(1);

  double timestepSize = 1.0;

  if (context.isNamed("M1SM")) {

    auto surfaceMeshID = interface.getMeshID("Surface_Mesh1");

    interface.setMeshVertices(surfaceMeshID, 2, surfaceCoords.data(), surface1VertexIDs.data());

    auto activationMeshID = interface.getMeshID("Activation_M1SM_Mesh");
    interface.setMeshVertices(activationMeshID, 1, signalCoords.data(), activationVertexIDs.data());

    auto stretchMeshID = interface.getMeshID("Stretch_M1SM_Mesh");
    interface.setMeshVertices(stretchMeshID, 1, signalCoords.data(), stretchVertexIDs.data());

    activationDataID    = interface.getDataID("Activation1", activationMeshID);
    stretchDataID       = interface.getDataID("Stretch1", stretchMeshID);
    displacement1DataID = interface.getDataID("Displacement1", surfaceMeshID);
    traction1DataID     = interface.getDataID("Traction1", surfaceMeshID);

  } else if (context.isNamed("M2SM")) {

    auto surfaceMeshID = interface.getMeshID("Surface_Mesh2");
    interface.setMeshVertices(surfaceMeshID, 2, surfaceCoords.data(), surface2VertexIDs.data());

    auto stretchMeshID = interface.getMeshID("Stretch_M2SM_Mesh");
    interface.setMeshVertices(stretchMeshID, 1, signalCoords.data(), stretchVertexIDs.data());

    stretchDataID       = interface.getDataID("Stretch2", stretchMeshID);
    displacement2DataID = interface.getDataID("Displacement2", surfaceMeshID);
    traction2DataID     = interface.getDataID("Traction2", surfaceMeshID);

  } else if (context.isNamed("M1")) {

    auto stretchMeshID = interface.getMeshID("Stretch_M1_Mesh");
    interface.setMeshVertices(stretchMeshID, 1, signalCoords.data(), stretchVertexIDs.data());

    auto crossStretchMeshID = interface.getMeshID("Stretch_M1_Cross_Mesh");
    interface.setMeshVertices(crossStretchMeshID, 1, signalCoords.data(), crossStretchVertexIDs.data());

    auto activationMeshID = interface.getMeshID("Activation_M1_Mesh");
    interface.setMeshVertices(activationMeshID, 1, signalCoords.data(), activationVertexIDs.data());

    stretchDataID      = interface.getDataID("Stretch1", stretchMeshID);
    crossStretchDataID = interface.getDataID("Stretch2", crossStretchMeshID);
    activationDataID   = interface.getDataID("Activation1", activationMeshID);

  } else {
    BOOST_TEST(context.isNamed("Tendon"));

    auto surface1MeshID = interface.getMeshID("SurfaceTendon_Mesh1");
    interface.setMeshVertices(surface1MeshID, 2, surfaceCoords.data(), surface1VertexIDs.data());

    auto surface2MeshID = interface.getMeshID("SurfaceTendon_Mesh2");
    interface.setMeshVertices(surface2MeshID, 2, surfaceCoords.data(), surface2VertexIDs.data());

    displacement1DataID = interface.getDataID("Displacement1", surface1MeshID);
    traction1DataID     = interface.getDataID("Traction1", surface1MeshID);

    displacement2DataID = interface.getDataID("Displacement2", surface2MeshID);
    traction2DataID     = interface.getDataID("Traction2", surface2MeshID);
  }

  std::cout << "Before initialize" << std::endl;
  interface.initialize();

  for (int timestep = 0; timestep < 2; ++timestep) {

    std::vector<double> tractions1{1.2, 3.4};
    std::vector<double> displacements1{4.2, 1.4};
    std::vector<double> tractions2{1.2, 3.7};
    std::vector<double> displacements2{4.1, 1.4};

    if (context.isNamed("M1SM")) {

      interface.writeBlockScalarData(displacement1DataID, 2, surface1VertexIDs.data(), displacements1.data());

    } else if (context.isNamed("Tendon")) {

      //interface.markActionFulfilled(precice::constants::actionWriteInitialData());

      std::vector<double> receivedDisplacements{0.0, 0.0};
      interface.readBlockScalarData(displacement1DataID, 2, surface1VertexIDs.data(), receivedDisplacements.data());
      BOOST_TEST(receivedDisplacements == displacements1, boost::test_tools::per_element());

    } else if (context.isNamed("M2SM")) {

    } else {
      BOOST_TEST(context.isNamed("M1"));
    }

    std::cout << "Before advance" << std::endl;
    interface.advance(timestepSize);
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Compositional

#endif // PRECICE_NO_MPI
