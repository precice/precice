#include "mesh/Vertex.hpp"
#include "testing/Testing.hpp"
#include "mesh/Edge.hpp"
#include <Eigen/Core>

using namespace precice::mesh;

BOOST_AUTO_TEST_SUITE(MeshTests)

BOOST_AUTO_TEST_CASE(Edges)
{
   Vertex v1 ( Eigen::Vector3d::Constant(0.0), 0 );
   Vertex v2 ( Eigen::Vector3d::Constant(1.0), 1 );
   Edge edge ( v1, v2, 0 );

   Eigen::VectorXd coords1 = edge.vertex(0).getCoords();
   Eigen::VectorXd coords2 = edge.vertex(1).getCoords();
   BOOST_TEST ( coords1 == Eigen::Vector3d::Constant(0.0) );
   BOOST_TEST ( coords2 == Eigen::Vector3d::Constant(1.0) );
}

BOOST_AUTO_TEST_SUITE_END() // Mesh
