#include <Eigen/Core>
#include <algorithm>
#include <math.h>
#include <ostream>
#include <vector>
#include "logging/Logger.hpp"
#include "math/differences.hpp"
#include "math/geometry.hpp"
#include "mesh/Data.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/RTree.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/impl/RTreeAdapter.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

using namespace precice;
using namespace precice::mesh;

namespace bg  = boost::geometry;
namespace bgi = boost::geometry::index;

BOOST_AUTO_TEST_SUITE(MeshTests)
BOOST_AUTO_TEST_SUITE(RTree)
BOOST_AUTO_TEST_SUITE(BGAdapters)

BOOST_AUTO_TEST_CASE(VectorAdapter)
{
  PRECICE_TEST(1_rank);
  Eigen::VectorXd vec = Eigen::Vector2d(1, 2);
  BOOST_TEST(bg::get<0>(vec) == 1);
  BOOST_TEST(bg::get<1>(vec) == 2);
  BOOST_TEST(bg::get<2>(vec) == 0);
  bg::set<1>(vec, 5);
  BOOST_TEST(bg::get<1>(vec) == 5);
}

BOOST_AUTO_TEST_CASE(VertexAdapter)
{
  PRECICE_TEST(1_rank);
  precice::mesh::Mesh mesh("MyMesh", 2, false, precice::testing::nextMeshID());
  auto &              v = mesh.createVertex(Eigen::Vector2d(1, 2));
  BOOST_TEST(bg::get<0>(v) == 1);
  BOOST_TEST(bg::get<1>(v) == 2);
  BOOST_TEST(bg::get<2>(v) == 0);
  bg::set<1>(v, 5);
  BOOST_TEST(bg::get<1>(v) == 5);
}

BOOST_AUTO_TEST_CASE(EdgeAdapter)
{
  PRECICE_TEST(1_rank);
  precice::mesh::Mesh mesh("MyMesh", 2, false, precice::testing::nextMeshID());
  auto &              v1 = mesh.createVertex(Eigen::Vector2d(1, 2));
  auto &              v2 = mesh.createVertex(Eigen::Vector2d(3, 4));
  auto &              e  = mesh.createEdge(v1, v2);
  BOOST_TEST((bg::get<0, 0>(e)) == 1.0);
  BOOST_TEST((bg::get<0, 1>(e)) == 2.0);
  BOOST_TEST((bg::get<0, 2>(e)) == 0.0);
  BOOST_TEST((bg::get<1, 0>(e)) == 3.0);
  BOOST_TEST((bg::get<1, 1>(e)) == 4.0);
  BOOST_TEST((bg::get<1, 2>(e)) == 0.0);
  bg::set<1, 1>(e, 5.0);
  BOOST_TEST((bg::get<1, 1>(e)) == 5.0);
}

BOOST_AUTO_TEST_CASE(TriangleAdapter)
{
  PRECICE_TEST(1_rank);
  precice::mesh::Mesh mesh("MyMesh", 3, false, precice::testing::nextMeshID());
  auto &              v1 = mesh.createVertex(Eigen::Vector3d(0, 2, 0));
  auto &              v2 = mesh.createVertex(Eigen::Vector3d(2, 1, 0));
  auto &              v3 = mesh.createVertex(Eigen::Vector3d(1, 0, 0));
  auto &              e1 = mesh.createEdge(v1, v2);
  auto &              e2 = mesh.createEdge(v2, v3);
  auto &              e3 = mesh.createEdge(v3, v1);
  auto &              t  = mesh.createTriangle(e1, e2, e3);

  std::vector<Eigen::VectorXd> vertices(t.begin(), t.end());
  std::vector<Eigen::VectorXd> refs{v1.getCoords(), v2.getCoords(), v3.getCoords()};
  BOOST_TEST(vertices.size() == refs.size());
  BOOST_TEST((std::is_permutation(
      vertices.begin(), vertices.end(),
      refs.begin(),
      [](const Eigen::VectorXd &lhs, const Eigen::VectorXd &rhs) {
        return precice::math::equals(lhs, rhs);
      })));
}

BOOST_AUTO_TEST_CASE(DistanceTestFlatSingleTriangle)
{
  PRECICE_TEST(1_rank);
  precice::mesh::Mesh mesh("MyMesh", 3, false, testing::nextMeshID());
  auto &              v1 = mesh.createVertex(Eigen::Vector3d(0, 0, 0));
  auto &              v2 = mesh.createVertex(Eigen::Vector3d(0, 1, 0));
  auto &              v3 = mesh.createVertex(Eigen::Vector3d(1, 0, 0));
  auto &              v4 = mesh.createVertex(Eigen::Vector3d(1, 1, 0));
  auto &              v5 = mesh.createVertex(Eigen::Vector3d(0.2, 0.2, 0));
  auto &              e1 = mesh.createEdge(v1, v2);
  auto &              e2 = mesh.createEdge(v2, v3);
  auto &              e3 = mesh.createEdge(v3, v1);
  auto &              t  = mesh.createTriangle(e1, e2, e3);

  BOOST_TEST(bg::comparable_distance(v1, v2) > 0.5);
  BOOST_TEST(bg::comparable_distance(v1, v1) < 0.01);
  BOOST_TEST(bg::comparable_distance(e1, v2) < 0.01);
  BOOST_TEST(bg::comparable_distance(e1, v3) > 0.2);
  BOOST_TEST(bg::comparable_distance(t, v3) < 0.1);
  BOOST_TEST(bg::comparable_distance(t, v4) > 0.2);
  BOOST_TEST(bg::comparable_distance(t, v5) < 0.01);
}

BOOST_AUTO_TEST_CASE(DistanceTestFlatDoubleTriangle)
{
  PRECICE_TEST(1_rank);
  precice::mesh::Mesh mesh("MyMesh", 3, false, testing::nextMeshID());
  auto &              lv1 = mesh.createVertex(Eigen::Vector3d(-1, 1, 0.1));
  auto &              lv2 = mesh.createVertex(Eigen::Vector3d(0, -1, 0));
  auto &              lv3 = mesh.createVertex(Eigen::Vector3d(-2, 0, -0.1));
  auto &              le1 = mesh.createEdge(lv1, lv2);
  auto &              le2 = mesh.createEdge(lv2, lv3);
  auto &              le3 = mesh.createEdge(lv3, lv1);
  auto &              lt  = mesh.createTriangle(le1, le2, le3);

  auto &rv1 = mesh.createVertex(Eigen::Vector3d(0, 1, 0.1));
  auto &rv2 = mesh.createVertex(Eigen::Vector3d(2, 0, -0.1));
  auto &rv3 = mesh.createVertex(Eigen::Vector3d(1, -1, 0));
  auto &re1 = mesh.createEdge(rv1, rv2);
  auto &re2 = mesh.createEdge(rv2, rv3);
  auto &re3 = mesh.createEdge(rv3, rv1);
  auto &rt  = mesh.createTriangle(re1, re2, re3);

  auto &v1 = mesh.createVertex(Eigen::Vector3d(-2, 1, 0));
  auto &v2 = mesh.createVertex(Eigen::Vector3d(2, -1, 0));
  auto &v3 = mesh.createVertex(Eigen::Vector3d(0, 0, 0));

  auto lt_v1 = bg::comparable_distance(lt, v1);
  auto lt_v2 = bg::comparable_distance(lt, v2);

  auto rt_v1 = bg::comparable_distance(rt, v1);
  auto rt_v3 = bg::comparable_distance(rt, v3);

  BOOST_TEST(precice::testing::equals(lt_v1, 0.5));
  BOOST_TEST(lt_v2 > 0);
  BOOST_TEST(rt_v1 > 1);
  BOOST_TEST(rt_v3 > 0);
}

BOOST_AUTO_TEST_CASE(DistanceTestFlatDoubleTriangleInsideOutside)
{
  PRECICE_TEST(1_rank);
  precice::mesh::Mesh mesh("MyMesh", 3, false, testing::nextMeshID());
  auto &              a = mesh.createVertex(Eigen::Vector3d(0, 0, 0));
  auto &              b = mesh.createVertex(Eigen::Vector3d(1, 0, 0));
  auto &              c = mesh.createVertex(Eigen::Vector3d(1, 1, 0));
  auto &              d = mesh.createVertex(Eigen::Vector3d(0, 1, 0));

  auto &ab = mesh.createEdge(a, b);
  auto &bd = mesh.createEdge(b, d);
  auto &da = mesh.createEdge(d, a);
  auto &dc = mesh.createEdge(d, c);
  auto &cb = mesh.createEdge(c, b);

  auto &lt = mesh.createTriangle(ab, bd, da);
  auto &rt = mesh.createTriangle(bd, dc, cb);
  BOOST_TEST_MESSAGE("Left  Triangle:" << lt);
  BOOST_TEST_MESSAGE("Right Triangle:" << rt);

  auto &lv = mesh.createVertex(Eigen::Vector3d(.25, .25, 0));
  auto &rv = mesh.createVertex(Eigen::Vector3d(.75, .75, 0));

  auto lv_lt = bg::comparable_distance(lv, lt);
  auto lv_rt = bg::comparable_distance(lv, rt);
  BOOST_TEST(precice::testing::equals(lv_lt, 0.0), lv_lt << " == 0.0");
  BOOST_TEST(!precice::testing::equals(lv_rt, 0.0), lv_rt << " != 0.0");

  auto rv_lt = bg::comparable_distance(rv, lt);
  auto rv_rt = bg::comparable_distance(rv, rt);
  BOOST_TEST(!precice::testing::equals(rv_lt, 0.0), rv_lt << " != 0.0");
  BOOST_TEST(precice::testing::equals(rv_rt, 0.0), rv_rt << " == 0.0");
}

BOOST_AUTO_TEST_CASE(DistanceTestSlopedTriangle)
{
  PRECICE_TEST(1_rank);
  precice::mesh::Mesh mesh("MyMesh", 3, false, testing::nextMeshID());
  auto &              v1 = mesh.createVertex(Eigen::Vector3d(0, 1, 0));
  auto &              v2 = mesh.createVertex(Eigen::Vector3d(1, 1, 1));
  auto &              v3 = mesh.createVertex(Eigen::Vector3d(0, 0, 1));
  auto &              v4 = mesh.createVertex(Eigen::Vector3d(0, 1, 1));
  auto &              v5 = mesh.createVertex(Eigen::Vector3d(1, 0, 0));
  auto &              e1 = mesh.createEdge(v1, v2);
  auto &              e2 = mesh.createEdge(v2, v3);
  auto &              e3 = mesh.createEdge(v3, v1);
  auto &              t  = mesh.createTriangle(e1, e2, e3);

  auto t_v4 = bg::comparable_distance(t, v4);
  auto v4_t = bg::comparable_distance(v4, t);
  BOOST_TEST(t_v4 > 0.01);
  BOOST_TEST(v4_t > 0.01);
  BOOST_TEST(v4_t == t_v4);

  auto t_v5 = bg::comparable_distance(t, v5);
  auto v5_t = bg::comparable_distance(v5, t);
  BOOST_TEST(t_v5 > 0.01);
  BOOST_TEST(v5_t > 0.01);
  BOOST_TEST(v5_t == t_v5);

  BOOST_TEST(v4_t < t_v5);
}

BOOST_AUTO_TEST_CASE(EnvelopeTriangleClockWise)
{
  PRECICE_TEST(1_rank);
  using precice::testing::equals;
  precice::mesh::Mesh mesh("MyMesh", 3, false, testing::nextMeshID());
  auto &              v1  = mesh.createVertex(Eigen::Vector3d(0, 1, 0));
  auto &              v2  = mesh.createVertex(Eigen::Vector3d(1, 1, 1));
  auto &              v3  = mesh.createVertex(Eigen::Vector3d(0, 0, 1));
  auto &              e1  = mesh.createEdge(v1, v2);
  auto &              e2  = mesh.createEdge(v2, v3);
  auto &              e3  = mesh.createEdge(v3, v1);
  auto &              t   = mesh.createTriangle(e1, e2, e3);
  auto                box = bg::return_envelope<precice::mesh::RTreeBox>(t);
  BOOST_TEST(equals(box.min_corner(), Eigen::Vector3d{0, 0, 0}));
  BOOST_TEST(equals(box.max_corner(), Eigen::Vector3d{1, 1, 1}));
}

BOOST_AUTO_TEST_CASE(EnvelopeTriangleCounterclockWise)
{
  PRECICE_TEST(1_rank);
  using precice::testing::equals;
  precice::mesh::Mesh mesh("MyMesh", 3, false, testing::nextMeshID());
  auto &              v1  = mesh.createVertex(Eigen::Vector3d(0, 1, 0));
  auto &              v2  = mesh.createVertex(Eigen::Vector3d(1, 1, 1));
  auto &              v3  = mesh.createVertex(Eigen::Vector3d(0, 0, 1));
  auto &              e1  = mesh.createEdge(v1, v3);
  auto &              e2  = mesh.createEdge(v3, v2);
  auto &              e3  = mesh.createEdge(v2, v1);
  auto &              t   = mesh.createTriangle(e1, e2, e3);
  auto                box = bg::return_envelope<precice::mesh::RTreeBox>(t);
  BOOST_TEST(equals(box.min_corner(), Eigen::Vector3d{0, 0, 0}));
  BOOST_TEST(equals(box.max_corner(), Eigen::Vector3d{1, 1, 1}));
}

BOOST_AUTO_TEST_SUITE_END() // BG Adapters
BOOST_AUTO_TEST_SUITE_END() // RTree
BOOST_AUTO_TEST_SUITE_END() // Mesh
