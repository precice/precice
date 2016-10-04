#include "VertexTest.hpp"
#include "mesh/Vertex.hpp"
#include "utils/Dimensions.hpp"
#include "utils/Parallel.hpp"
#include "math/math.hpp"
#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::mesh::tests::VertexTest)

namespace precice {
namespace mesh {
namespace tests {

logging::Logger VertexTest:: _log ( "precice::mesh::VertexTest" );

VertexTest:: VertexTest ()
:
  TestCase ("mesh::VertexTest")
{}

void VertexTest:: run ()
{
  PRECICE_MASTER_ONLY {
    testMethod ( test );
  }
}

void VertexTest:: test ()
{
  preciceTrace ( "test()" );
  
  Vertex vertex ( Eigen::Vector3d::Constant(1.0), 0 );

  Eigen::Vector3d coords = vertex.getCoords();
  validate ( math::equals(coords, Eigen::Vector3d::Constant(1.0)) );

  int id = vertex.getID ();
  validateEquals ( id, 0 );

  Eigen::Vector3d normal = vertex.getNormal();
  validate ( math::equals(normal, Eigen::Vector3d::Zero()) );

  void* mesh = static_cast<void*> ( vertex.mesh() );
  // Can be replaced by nullptr as soon as we have C++11 available.
  validateEquals ( mesh, static_cast<void*>(nullptr) );
}

}}} // namespace precice, mesh, tests
