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
BOOST_AUTO_TEST_CASE(EdgeEquality)
{
   Vertex v1 (Eigen::Vector3d(0,0,0), 0);
   Vertex v2 (Eigen::Vector3d(0,0,1), 0);
   Vertex v3 (Eigen::Vector3d(0,0,2), 0);
   Edge edge1(v1, v2, 0);
   Edge edge2(v2, v1, 1);
   Edge edge3(v1, v3, 0);
   Edge edge4(v1, v3, 0);
   edge4.setNormal(Eigen::Vector3d(Eigen::Vector3d(0,1,0)));
   BOOST_TEST(edge1 == edge2);
   BOOST_TEST(edge1 != edge3);
   BOOST_TEST(edge3 != edge4);
}

BOOST_AUTO_TEST_SUITE_END() // Mesh
