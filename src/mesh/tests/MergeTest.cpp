#include "mesh/Edge.hpp"
#include "mesh/Group.hpp"
#include "mesh/Merge.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Quad.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Vertex.hpp"
#include "testing/Testing.hpp"

using namespace precice::mesh;

BOOST_AUTO_TEST_SUITE(MeshTests)

BOOST_AUTO_TEST_CASE(MergeTest)
{
  // Create visitables
  int                 dim = 3;
  precice::mesh::Mesh mesh("MyMesh", dim, false);
  Vertex &            v0 = mesh.createVertex(Eigen::Vector3d::Zero());
  Vertex &            v1 = mesh.createVertex(Eigen::Vector3d::Zero());
  Vertex &            v2 = mesh.createVertex(Eigen::Vector3d::Zero());
  Vertex &            v3 = mesh.createVertex(Eigen::Vector3d::Zero());
  Edge &              e0 = mesh.createEdge(v0, v1);
  Edge &              e1 = mesh.createEdge(v1, v2);
  Edge &              e2 = mesh.createEdge(v2, v0);
  Edge &              e3 = mesh.createEdge(v2, v3);
  Edge &              e4 = mesh.createEdge(v3, v0);
  Triangle &          t1 = mesh.createTriangle(e0, e1, e2);
  Triangle &          t2 = mesh.createTriangle(e2, e1, e0);
  Quad &              q0 = mesh.createQuad(e0, e1, e3, e4);

  Group group;
  group.add(v1);
  group.add(v2);
  group.add(v3);
  group.add(v2);
  group.add(v1);
  group.add(v3);
  group.add(v2);

  group.add(e1);
  group.add(e1);
  group.add(e2);
  group.add(e3);
  group.add(e3);
  group.add(e4);

  group.add(t1);
  group.add(t1);
  group.add(t2);
  group.add(t2);

  group.add(q0);
  group.add(q0);

  Merge merge;
  merge(group);
  BOOST_TEST(merge.content().vertices().size() == 3);
  BOOST_TEST(merge.content().edges().size() == 4);
  BOOST_TEST(merge.content().triangles().size() == 2);
  BOOST_TEST(merge.content().quads().size() == 1);
  BOOST_TEST(merge.content().size() == 10);
  merge(group); // Shouldn't change anything
  BOOST_TEST(merge.content().vertices().size() == 3);
  BOOST_TEST(merge.content().edges().size() == 4);
  BOOST_TEST(merge.content().triangles().size() == 2);
  BOOST_TEST(merge.content().quads().size() == 1);
  BOOST_TEST(merge.content().size() == 10);
}

BOOST_AUTO_TEST_SUITE_END() // Mesh
