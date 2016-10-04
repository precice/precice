#include "EdgeTest.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Edge.hpp"
#include "utils/Parallel.hpp"
#include "utils/Dimensions.hpp"

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::mesh::tests::EdgeTest)

namespace precice {
namespace mesh {
namespace tests {

logging::Logger EdgeTest:: _log ( "precice::mesh::EdgeTest" );

EdgeTest:: EdgeTest ()
:
   TestCase ("mesh::EdgeTest")
{}

void EdgeTest:: run ()
{
   PRECICE_MASTER_ONLY {
      testMethod ( test );
   }
}

void EdgeTest:: test ()
{
   preciceTrace ( "test()" );

   Vertex v1 ( Eigen::Vector3d::Constant(0.0), 0 );
   Vertex v2 ( Eigen::Vector3d::Constant(1.0), 1 );
//   VertexTuple vertices = { Vertex(Vector(0.0)), Vertex(Vector(1.0)) };
   Edge edge ( v1, v2, 0 );

   Eigen::VectorXd coords1 = edge.vertex(0).getCoords();
   Eigen::VectorXd coords2 = edge.vertex(1).getCoords();
   validate ( coords1 == Eigen::Vector3d::Constant(0.0) );
   validate ( coords2 == Eigen::Vector3d::Constant(1.0) );

   Edge edge2 ( v1, v2, 1 );
}

}}} // namespace precice, mesh, tests
