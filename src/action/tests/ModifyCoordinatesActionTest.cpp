#include "ModifyCoordinatesActionTest.hpp"
#include "action/ModifyCoordinatesAction.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Edge.hpp"
#include "utils/Parallel.hpp"
#include "math/math.hpp"

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::action::tests::ModifyCoordinatesActionTest)

namespace precice {
namespace action {
namespace tests {

logging::Logger ModifyCoordinatesActionTest::
_log ( "precice::action::tests::ModifyCoordinatesActionTest" );

ModifyCoordinatesActionTest:: ModifyCoordinatesActionTest ()
  :
  TestCase ( "precice::action::tests::ModifyCoordinatesActionTest" )
{}

void ModifyCoordinatesActionTest:: run ()
{
  PRECICE_MASTER_ONLY {
    testMethod ( testAddToCoordinates );
    testMethod ( testSubtractFromCoordinates );
  }
}

void ModifyCoordinatesActionTest:: testAddToCoordinates ()
{
  TRACE();
  using namespace mesh;
  PtrMesh mesh ( new Mesh("Mesh", 2, false) );
  PtrData data = mesh->createData ( "test-data", 2 );
  int dataID = data->getID ();
  Vertex& v0 = mesh->createVertex ( Eigen::Vector2d::Zero() );
  Vertex& v1 = mesh->createVertex ( Eigen::Vector2d::Constant(1.0) );
  Edge & edge = mesh->createEdge ( v0, v1 );
  mesh->computeState();
  mesh->allocateDataValues ();
  auto& values = data->values ();

  // Create ApplyDisplacementsMeshAction
  action::ModifyCoordinatesAction modifyCoordinates (
  action::ModifyCoordinatesAction::ALWAYS_PRIOR, dataID, mesh,
  action::ModifyCoordinatesAction::ADD_TO_COORDINATES_MODE );

  // Validate coordinates of mesh before modifying it
  validate ( math::equals(v0.getCoords(), Eigen::Vector2d::Constant(0.0)) );
  validate ( math::equals(v1.getCoords(), Eigen::Vector2d::Constant(1.0)) );
  Eigen::Vector2d normalizedNormal(0.5, -0.5);
  normalizedNormal.normalize();
  validate ( math::equals(edge.getNormal(), normalizedNormal) );

  // Set displacements
  values.segment(v0.getID()*2, 2) = Eigen::VectorXd::Constant(2, 2.0);
  values.segment(v1.getID()*2, 2) = Eigen::VectorXd::Constant(2, -2.0);
  
  // Apply displacements to  node coordinates
  modifyCoordinates.performAction(0.0, 0.0, 0.0, 0.0);

  // Validate coordinates of mesh after modifying it
  validate(math::equals(v0.getCoords(), Eigen::Vector2d::Constant(2.0)));
  validate(math::equals(v1.getCoords(), Eigen::Vector2d::Constant(-1.0)));
  normalizedNormal << -0.5, 0.5;
  normalizedNormal.normalize();
  validate(math::equals(edge.getNormal(), normalizedNormal));
}

void ModifyCoordinatesActionTest:: testSubtractFromCoordinates ()
{
  TRACE();
  using namespace mesh;
  // Create geometryContext by faking a geometry but not using it to create
  // the mesh. The mesh is created by hand, such that references to the vertices
  // to be displaced are obtained.
  PtrMesh mesh ( new Mesh("Mesh", 2, false) );
  PtrData data = mesh->createData ( "test-data", 2 );
  int dataID = data->getID ();
  Vertex& v0 = mesh->createVertex ( Eigen::Vector2d::Constant(0.0) );
  Vertex& v1 = mesh->createVertex ( Eigen::Vector2d::Constant(1.0) );
  Edge& edge = mesh->createEdge ( v0, v1 );
  mesh->computeState();
  mesh->allocateDataValues ();
  auto& values = data->values ();

  // Create ApplyDisplacementsMeshAction
  action::ModifyCoordinatesAction modifyCoordinates (
    action::ModifyCoordinatesAction::ALWAYS_PRIOR, dataID, mesh,
    action::ModifyCoordinatesAction::SUBTRACT_FROM_COORDINATES_MODE );
//   modifyCoordinates.loadMeshContext ( meshContext );

  // Validate coordinates of mesh before modifying it
  validate ( math::equals(v0.getCoords(), Eigen::Vector2d::Constant(0.0)) );
  validate ( math::equals(v1.getCoords(), Eigen::Vector2d::Constant(1.0)) );
  Eigen::Vector2d normalizedNormal(0.5, -0.5);
  normalizedNormal.normalize();
  validate ( math::equals(edge.getNormal(), normalizedNormal) );

  // Set displacements
  values.segment(v0.getID()*2, 2) = Eigen::VectorXd::Constant(2, -2.0);
  values.segment(v1.getID()*2, 2) = Eigen::VectorXd::Constant(2, 2.0);
  //tarch::la::slice<2>(values,v0.getID()*2) = Eigen::Vector2d(-2.0);
  //tarch::la::slice<2>(values,v1.getID()*2) = Eigen::Vector2d(2.0);

  // Apply displacements to  node coordinates
  modifyCoordinates.performAction(0.0, 0.0, 0.0, 0.0);

  // Validate coordinates of mesh after modifying it
  validate ( math::equals(v0.getCoords(), Eigen::Vector2d::Constant(2.0)) );
  validate ( math::equals(v1.getCoords(), Eigen::Vector2d::Constant(-1.0)) );
  normalizedNormal << -0.5, 0.5;
  normalizedNormal.normalize();
  validate ( math::equals(edge.getNormal(), normalizedNormal) );
}

}}} // namespace precice, action, tests
