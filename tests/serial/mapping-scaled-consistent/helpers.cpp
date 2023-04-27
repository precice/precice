#ifndef PRECICE_NO_MPI

#include "helpers.hpp"
#include "math/geometry.hpp"
#include "precice/SolverInterface.hpp"
#include "precice/impl/SolverInterfaceImpl.hpp"
#include "testing/Testing.hpp"

void testQuadMappingScaledConsistent(const std::string configFile, const TestContext &context)
{
  using Eigen::Vector3d;

  const double z = 0.3;

  // MeshOne
  Vector3d coordOneA{0.0, 0.0, z};
  Vector3d coordOneB{1.0, 0.0, z};
  Vector3d coordOneC{0.999999999, 1.0, z}; // Forces diagonal 0-2 to be shorter.
  Vector3d coordOneD{0.0, 1.0, z};
  double   valOneA = 1.0;
  double   valOneB = 3.0;
  double   valOneC = 5.0;
  double   valOneD = 7.0;

  // MeshTwo
  Vector3d coordTwoA{0.0, 0.0, z + 0.1};               // Maps to vertex A
  Vector3d coordTwoB{0.0, 0.5, z - 0.01};              // Maps to edge AD
  Vector3d coordTwoC{2.0 / 3.0, 1.0 / 3.0, z + 0.001}; // Maps to triangle ABC

  double expectedIntegral = precice::math::geometry::triangleArea(coordOneA, coordOneB, coordOneC) * (valOneA + valOneB + valOneC) / 3.0 +
                            precice::math::geometry::triangleArea(coordOneA, coordOneC, coordOneD) * (valOneA + valOneC + valOneD) / 3.0;

  if (context.isNamed("SolverOne")) {
    precice::SolverInterface interface("SolverOne", configFile, 0, 1);
    // namespace is required because we are outside the fixture
    auto meshOneID = "MeshOne";

    // Setup mesh one.
    int idA = interface.setMeshVertex(meshOneID, coordOneA.data());
    int idB = interface.setMeshVertex(meshOneID, coordOneB.data());
    int idC = interface.setMeshVertex(meshOneID, coordOneC.data());
    int idD = interface.setMeshVertex(meshOneID, coordOneD.data());

    interface.setMeshQuad(meshOneID, idA, idB, idC, idD);

    auto &mesh = testing::WhiteboxAccessor::impl(interface).mesh("MeshOne");
    BOOST_REQUIRE(mesh.vertices().size() == 4);
    BOOST_REQUIRE(mesh.edges().empty());
    BOOST_REQUIRE(mesh.triangles().size() == 2);

    // Initialize, thus sending the mesh.
    interface.initialize();
    double maxDt = interface.getMaxTimeStepSize();
    BOOST_TEST(mesh.edges().size() == 5);
    BOOST_TEST(mesh.triangles().size() == 2);
    BOOST_TEST(interface.isCouplingOngoing(), "Sending participant should have to advance once!");

    // Write the data to be send.
    auto   dataAID  = "DataOne";
    int    ids[]    = {idA, idB, idC, idD};
    double values[] = {valOneA, valOneB, valOneC, valOneD};
    interface.writeData(meshOneID, dataAID, ids, values);

    // Advance, thus send the data to the receiving partner.
    interface.advance(maxDt);
    BOOST_TEST(!interface.isCouplingOngoing(), "Sending participant should have to advance once!");
    interface.finalize();
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    precice::SolverInterface interface("SolverTwo", configFile, 0, 1);
    // namespace is required because we are outside the fixture
    auto meshTwoID = "MeshTwo";

    // Setup receiving mesh.
    int idA = interface.setMeshVertex(meshTwoID, coordTwoA.data());
    int idB = interface.setMeshVertex(meshTwoID, coordTwoB.data());
    int idC = interface.setMeshVertex(meshTwoID, coordTwoC.data());

    interface.setMeshTriangle(meshTwoID, idA, idB, idC);

    // Initialize, thus receive the data and map.
    interface.initialize();
    double maxDt = interface.getMaxTimeStepSize();
    BOOST_TEST(interface.isCouplingOngoing(), "Receiving participant should have to advance once!");

    // Read the mapped data from the mesh.
    auto   dataAID = "DataOne";
    int    ids[]   = {idA, idB, idC};
    double values[3];
    interface.readData(meshTwoID, dataAID, ids, maxDt, values);

    double calculatedIntegral = precice::math::geometry::triangleArea(coordTwoA, coordTwoB, coordTwoC) * (values[0] + values[1] + values[2]) / 3.0;
    BOOST_TEST(expectedIntegral == calculatedIntegral);

    // Verify that there is only one time step necessary.
    interface.advance(maxDt);
    BOOST_TEST(!interface.isCouplingOngoing(), "Receiving participant should have to advance once!");
    interface.finalize();
  }
}

void testQuadMappingScaledConsistentVolumetric(const std::string configFile, const TestContext &context)
{
  using Eigen::Vector2d;

  // MeshOne
  Vector2d coordOneA{0.0, 0.0};
  Vector2d coordOneB{1.0, 0.0};
  Vector2d coordOneC{1.0, 1.0};
  Vector2d coordOneD{0.0, 1.0};
  Vector2d coordOneExtra{0.5, 1.0};
  double   valOneA     = 4.0;
  double   valOneB     = 3.0;
  double   valOneC     = 2.0;
  double   valOneD     = 1.0;
  double   valOneExtra = 1.0;

  // MeshTwo
  Vector2d coordTwoA{0.0, 0.0};
  Vector2d coordTwoB{1.0, 0.0};
  Vector2d coordTwoC{1.0, 1.0};
  Vector2d coordTwoD{0.0, 1.0};

  double expectedIntegral = 7.0 / 3;

  if (context.isNamed("SolverOne")) {
    precice::SolverInterface interface("SolverOne", configFile, 0, 1);
    // namespace is required because we are outside the fixture
    auto meshOneID = "MeshOne";

    // Setup mesh one.
    int idA     = interface.setMeshVertex(meshOneID, coordOneA.data());
    int idB     = interface.setMeshVertex(meshOneID, coordOneB.data());
    int idC     = interface.setMeshVertex(meshOneID, coordOneC.data());
    int idD     = interface.setMeshVertex(meshOneID, coordOneD.data());
    int idExtra = interface.setMeshVertex(meshOneID, coordOneExtra.data());

    interface.setMeshTriangle(meshOneID, idA, idB, idExtra);
    interface.setMeshTriangle(meshOneID, idA, idD, idExtra);
    interface.setMeshTriangle(meshOneID, idB, idC, idExtra);

    auto &mesh = testing::WhiteboxAccessor::impl(interface).mesh("MeshOne");
    BOOST_REQUIRE(mesh.vertices().size() == 5);
    BOOST_REQUIRE(mesh.triangles().size() == 3);

    // Initialize, thus sending the mesh.
    interface.initialize();
    double maxDt = interface.getMaxTimeStepSize();
    BOOST_TEST(interface.isCouplingOngoing(), "Sending participant should have to advance once!");

    // Write the data to be send.
    auto   dataAID  = "DataOne";
    int    ids[]    = {idA, idB, idC, idD, idExtra};
    double values[] = {valOneA, valOneB, valOneC, valOneD, valOneExtra};
    interface.writeData(meshOneID, dataAID, ids, values);

    // Advance, thus send the data to the receiving partner.
    interface.advance(maxDt);
    BOOST_TEST(!interface.isCouplingOngoing(), "Sending participant should have to advance once!");
    interface.finalize();
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    precice::SolverInterface interface("SolverTwo", configFile, 0, 1);
    // namespace is required because we are outside the fixture
    auto meshTwoID = "MeshTwo";

    // Setup receiving mesh.
    int idA = interface.setMeshVertex(meshTwoID, coordTwoA.data());
    int idB = interface.setMeshVertex(meshTwoID, coordTwoB.data());
    int idC = interface.setMeshVertex(meshTwoID, coordTwoC.data());
    int idD = interface.setMeshVertex(meshTwoID, coordTwoD.data());

    interface.setMeshEdge(meshTwoID, idA, idB);
    interface.setMeshEdge(meshTwoID, idB, idC);
    interface.setMeshEdge(meshTwoID, idA, idC);

    interface.setMeshTriangle(meshTwoID, idA, idB, idC);
    interface.setMeshTriangle(meshTwoID, idA, idD, idC);

    // Initialize, thus receive the data and map.
    interface.initialize();
    double maxDt = interface.getMaxTimeStepSize();
    BOOST_TEST(interface.isCouplingOngoing(), "Receiving participant should have to advance once!");

    // Read the mapped data from the mesh.
    auto dataAID = "DataOne";

    int    ids[] = {idA, idB, idC, idD};
    double values[4];
    interface.readData(meshTwoID, dataAID, ids, maxDt, values);

    double calculatedIntegral = precice::math::geometry::triangleArea(coordTwoA, coordTwoB, coordTwoC) * (values[0] + values[1] + values[2]) / 3.0 +
                                precice::math::geometry::triangleArea(coordTwoA, coordTwoD, coordTwoC) * (values[0] + values[1] + values[2]) / 3.0;
    BOOST_TEST(expectedIntegral == calculatedIntegral);
    BOOST_TEST(values[0] = valOneA * 8.0 / 7);

    // Verify that there is only one time step necessary.
    interface.advance(maxDt);
    BOOST_TEST(!interface.isCouplingOngoing(), "Receiving participant should have to advance once!");
    interface.finalize();
  }
}

void testTetraScaledConsistentVolumetric(const std::string configFile, const TestContext &context)
{
  using Eigen::Vector3d;

  // MeshOne
  Vector3d coordOneA{0.0, 0.0, 0.0};
  Vector3d coordOneB{0.5, 1.0, 0.0};
  Vector3d coordOneC{1.0, 0.0, 0.0};
  Vector3d coordOneD{0.5, 0.0, 1.0};
  Vector3d coordOneExtra{0.5, 0.0, 0.0};
  double   valOneA     = 1.0;
  double   valOneB     = 2.0;
  double   valOneC     = 3.0;
  double   valOneD     = 4.0;
  double   valOneExtra = 5.0;

  // MeshTwo
  Vector3d coordTwoA{0.0, 0.0, 0.0};
  Vector3d coordTwoB{0.5, 1.0, 0.0};
  Vector3d coordTwoC{1.0, 0.0, 0.0};
  Vector3d coordTwoD{0.5, 0.0, 1.0};

  double expectedIntegral = 6.5 / 12;

  if (context.isNamed("SolverOne")) {
    precice::SolverInterface interface("SolverOne", configFile, 0, 1);
    // namespace is required because we are outside the fixture
    auto meshOneID = "MeshOne";

    // Setup mesh one.
    int idA     = interface.setMeshVertex(meshOneID, coordOneA.data());
    int idB     = interface.setMeshVertex(meshOneID, coordOneB.data());
    int idC     = interface.setMeshVertex(meshOneID, coordOneC.data());
    int idD     = interface.setMeshVertex(meshOneID, coordOneD.data());
    int idExtra = interface.setMeshVertex(meshOneID, coordOneExtra.data());

    interface.setMeshTetrahedron(meshOneID, idA, idB, idD, idExtra);
    interface.setMeshTetrahedron(meshOneID, idC, idB, idD, idExtra);

    auto &mesh = testing::WhiteboxAccessor::impl(interface).mesh("MeshOne");
    BOOST_REQUIRE(mesh.vertices().size() == 5);
    BOOST_REQUIRE(mesh.tetrahedra().size() == 2);

    // Initialize, thus sending the mesh.
    interface.initialize();
    double maxDt = interface.getMaxTimeStepSize();
    BOOST_TEST(interface.isCouplingOngoing(), "Sending participant should have to advance once!");

    // Write the data to be send.
    auto   dataAID  = "DataOne";
    int    ids[]    = {idA, idB, idC, idD, idExtra};
    double values[] = {valOneA, valOneB, valOneC, valOneD, valOneExtra};
    interface.writeData(meshOneID, dataAID, ids, values);

    // Advance, thus send the data to the receiving partner.
    interface.advance(maxDt);
    BOOST_TEST(!interface.isCouplingOngoing(), "Sending participant should have to advance once!");
    interface.finalize();
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    precice::SolverInterface interface("SolverTwo", configFile, 0, 1);
    // namespace is required because we are outside the fixture
    auto meshTwoID = "MeshTwo";

    // Setup receiving mesh.
    int idA = interface.setMeshVertex(meshTwoID, coordTwoA.data());
    int idB = interface.setMeshVertex(meshTwoID, coordTwoB.data());
    int idC = interface.setMeshVertex(meshTwoID, coordTwoC.data());
    int idD = interface.setMeshVertex(meshTwoID, coordTwoD.data());

    interface.setMeshTetrahedron(meshTwoID, idA, idB, idC, idD);

    // Initialize, thus receive the data and map.
    interface.initialize();
    double maxDt = interface.getMaxTimeStepSize();
    BOOST_TEST(interface.isCouplingOngoing(), "Receiving participant should have to advance once!");

    // Read the mapped data from the mesh.
    auto   dataAID = "DataOne";
    double values[4];
    int    ids[] = {idA, idB, idC, idD};
    interface.readData(meshTwoID, dataAID, ids, maxDt, values);

    double calculatedIntegral = precice::math::geometry::tetraVolume(coordTwoA, coordTwoB, coordTwoC, coordTwoD) * (values[0] + values[1] + values[2] + values[3]) / 4.0;
    BOOST_TEST(expectedIntegral == calculatedIntegral);
    BOOST_TEST(values[0] == valOneA * 1.3);

    // Verify that there is only one time step necessary.
    interface.advance(maxDt);
    BOOST_TEST(!interface.isCouplingOngoing(), "Receiving participant should have to advance once!");
    interface.finalize();
  }
}

#endif
