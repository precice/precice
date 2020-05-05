#include "math/geometry.hpp"
#include "mesh/RTree.hpp"
#include "mesh/impl/RTree.hpp"
#include "mesh/impl/RTreeAdapter.hpp"
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

BOOST_AUTO_TEST_CASE(QuadAdapter)
{
  PRECICE_TEST(1_rank);
  precice::mesh::Mesh mesh("MyMesh", 3, false, testing::nextMeshID());
  auto &              v1 = mesh.createVertex(Eigen::Vector3d(0, 2, 0));
  auto &              v2 = mesh.createVertex(Eigen::Vector3d(2, 1, 0));
  auto &              v3 = mesh.createVertex(Eigen::Vector3d(3, 0, 0));
  auto &              v4 = mesh.createVertex(Eigen::Vector3d(1, 0, 0));
  auto &              e1 = mesh.createEdge(v1, v2);
  auto &              e2 = mesh.createEdge(v2, v3);
  auto &              e3 = mesh.createEdge(v3, v4);
  auto &              e4 = mesh.createEdge(v4, v1);
  auto &              t  = mesh.createQuad(e1, e2, e3, e4);

  std::vector<Eigen::VectorXd> vertices(t.begin(), t.end());
  std::vector<Eigen::VectorXd> refs{v1.getCoords(), v2.getCoords(), v3.getCoords(), v4.getCoords()};
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

struct MeshFixture {
  MeshFixture()
      : mesh("MyMesh", 3, false, testing::nextMeshID())
  {
    auto &v1 = mesh.createVertex(Eigen::Vector3d(0, 2, 0));
    auto &v2 = mesh.createVertex(Eigen::Vector3d(2, 1, 0));
    auto &v3 = mesh.createVertex(Eigen::Vector3d(3, 0, 0));
    auto &v4 = mesh.createVertex(Eigen::Vector3d(1, 0, 0));
    // Quad Borders
    auto &e1 = mesh.createEdge(v1, v2);
    auto &e2 = mesh.createEdge(v2, v3);
    auto &e3 = mesh.createEdge(v3, v4);
    auto &e4 = mesh.createEdge(v4, v1);
    // Diagonal
    auto &e5 = mesh.createEdge(v2, v4);
    // Triangles
    mesh.createTriangle(e1, e5, e4);
    mesh.createTriangle(e2, e3, e5);
    // Quad
    mesh.createQuad(e1, e2, e3, e4);

    // Check the Mesh
    BOOST_TEST(mesh.vertices().size() == 4);
    BOOST_TEST(mesh.edges().size() == 5);
    BOOST_TEST(mesh.triangles().size() == 2);
    BOOST_TEST(mesh.quads().size() == 1);
  }

  Mesh mesh;

  const int vertex_cnt    = 4;
  const int edge_cnt      = 5;
  const int triangle_cnt  = 2;
  const int quad_cnt      = 1;
  const int primitive_cnt = 4 + 5 + 2 + 1;
};

BOOST_AUTO_TEST_CASE(Query_2D)
{
  PRECICE_TEST(1_rank);
  PtrMesh mesh(new precice::mesh::Mesh("MyMesh", 2, false, testing::nextMeshID()));
  mesh->createVertex(Eigen::Vector2d(0, 0));
  mesh->createVertex(Eigen::Vector2d(0, 1));
  auto &v1 = mesh->createVertex(Eigen::Vector2d(1, 0));
  auto &v2 = mesh->createVertex(Eigen::Vector2d(1, 1));
  mesh->createEdge(v1, v2);

  {
    auto tree = rtree::getVertexRTree(mesh);

    BOOST_TEST(tree->size() == 4);

    Eigen::VectorXd     searchVector(Eigen::Vector2d(0.2, 0.8));
    std::vector<size_t> results;

    tree->query(bgi::nearest(searchVector, 1), std::back_inserter(results));

    BOOST_TEST(results.size() == 1);
    BOOST_TEST(mesh->vertices()[results[0]].getCoords() == Eigen::Vector2d(0, 1));
  }

  {
    auto tree = rtree::getPrimitiveRTree(mesh);

    BOOST_TEST(tree->size() == 5);

    Eigen::VectorXd                         searchVector(Eigen::Vector2d(0.2, 0.8));
    std::vector<PrimitiveRTree::value_type> results;

    tree->query(bgi::nearest(searchVector, 1), std::back_inserter(results));

    BOOST_TEST(results.size() == 1);
    auto pi = results.front().second;
    BOOST_TEST(pi.type == Primitive::Vertex);
    BOOST_TEST(pi.index < mesh->vertices().size());
    BOOST_TEST(mesh->vertices()[pi.index].getCoords() == Eigen::Vector2d(0, 1));
  }
}

BOOST_AUTO_TEST_CASE(Query_3D)
{
  PRECICE_TEST(1_rank);
  PtrMesh mesh(new precice::mesh::Mesh("MyMesh", 3, false, precice::testing::nextMeshID()));
  mesh->createVertex(Eigen::Vector3d(0, 0, 0));
  mesh->createVertex(Eigen::Vector3d(0, 0, 1));
  mesh->createVertex(Eigen::Vector3d(0, 1, 0));
  mesh->createVertex(Eigen::Vector3d(0, 1, 1));
  mesh->createVertex(Eigen::Vector3d(1, 0, 0));
  mesh->createVertex(Eigen::Vector3d(1, 0, 1));
  auto &v1 = mesh->createVertex(Eigen::Vector3d(1, 1, 0));
  auto &v2 = mesh->createVertex(Eigen::Vector3d(1, 1, 1));
  mesh->createEdge(v1, v2);

  {
    auto tree = rtree::getVertexRTree(mesh);

    BOOST_TEST(tree->size() == 8);

    Eigen::VectorXd     searchVector(Eigen::Vector3d(0.8, 0.0, 0.8));
    std::vector<size_t> results;

    tree->query(bgi::nearest(searchVector, 1), std::back_inserter(results));

    BOOST_TEST(results.size() == 1);
    BOOST_TEST(mesh->vertices()[results[0]].getCoords() == Eigen::Vector3d(1, 0, 1));
  }

  {
    auto tree = rtree::getPrimitiveRTree(mesh);

    BOOST_TEST(tree->size() == 9);

    Eigen::VectorXd                         searchVector(Eigen::Vector3d(1.8, 0.0, 0.8));
    std::vector<PrimitiveRTree::value_type> results;

    tree->query(bgi::nearest(searchVector, 1), std::back_inserter(results));

    BOOST_TEST(results.size() == 1);
    auto pi = results.front().second;
    BOOST_TEST(pi.type == Primitive::Vertex);
    BOOST_TEST(pi.index < mesh->vertices().size());
    BOOST_TEST(mesh->vertices()[pi.index].getCoords() == Eigen::Vector3d(1, 0, 1));
  }

  {
    auto tree = rtree::getPrimitiveRTree(mesh);

    BOOST_TEST(tree->size() == 9);

    Eigen::VectorXd searchVector((v1.getCoords() + v2.getCoords()) / 2);
    searchVector += Eigen::Vector3d(0.001, -0.03, 0.005); // "noise"
    std::vector<PrimitiveRTree::value_type> results;

    tree->query(bgi::nearest(searchVector, 1), std::back_inserter(results));

    BOOST_TEST(results.size() == 1);
    auto pi = results.front().second;
    BOOST_TEST(pi.type == Primitive::Edge);
    BOOST_TEST(pi.index == 0);
  }
}

BOOST_AUTO_TEST_CASE(Query_3D_FullMesh)
{
  PRECICE_TEST(1_rank);
  PtrMesh      mesh(new precice::mesh::Mesh("MyMesh", 3, false, precice::testing::nextMeshID()));
  const double z1  = 0.1;
  const double z2  = -0.1;
  auto &       v00 = mesh->createVertex(Eigen::Vector3d(0, 0, 0));
  auto &       v01 = mesh->createVertex(Eigen::Vector3d(0, 1, 0));
  auto &       v10 = mesh->createVertex(Eigen::Vector3d(1, 0, z1));
  auto &       v11 = mesh->createVertex(Eigen::Vector3d(1, 1, z1));
  auto &       v20 = mesh->createVertex(Eigen::Vector3d(2, 0, z2));
  auto &       v21 = mesh->createVertex(Eigen::Vector3d(2, 1, z2));
  auto &       ell = mesh->createEdge(v00, v01);
  auto &       elt = mesh->createEdge(v01, v11);
  auto &       elr = mesh->createEdge(v11, v10);
  auto &       elb = mesh->createEdge(v10, v00);
  auto &       eld = mesh->createEdge(v00, v11);
  auto &       erl = elr;
  auto &       ert = mesh->createEdge(v11, v21);
  auto &       err = mesh->createEdge(v21, v20);
  auto &       erb = mesh->createEdge(v20, v10);
  auto &       erd = mesh->createEdge(v10, v21);
  auto &       tlt = mesh->createTriangle(ell, elt, eld);
  auto &       tlb = mesh->createTriangle(eld, elb, elr);
  auto &       trt = mesh->createTriangle(erl, ert, erd);
  auto &       trb = mesh->createTriangle(erd, erb, err);

  {
    auto tree = rtree::getVertexRTree(mesh);

    BOOST_TEST(tree->size() == 6);

    Eigen::VectorXd     searchVector(Eigen::Vector3d(0.8, 0.0, 0.8));
    std::vector<size_t> results;

    tree->query(bgi::nearest(searchVector, 1), std::back_inserter(results));

    BOOST_TEST_INFO(results);
    BOOST_TEST(results.size() == 1);
  }
  {
    auto tree = rtree::getEdgeRTree(mesh);

    BOOST_TEST(tree->size() == 9);

    Eigen::VectorXd  searchVector(Eigen::Vector3d(0.8, 0.5, 0.0));
    std::set<size_t> results;

    tree->query(bgi::nearest(searchVector, 2), std::inserter(results, results.begin()));

    BOOST_TEST_INFO(results);
    BOOST_TEST(results.size() == 2);
    BOOST_TEST(results.count(eld.getID()) == 1);
    BOOST_TEST(results.count(elr.getID()) == 1);
  }
  {
    auto tree = rtree::getTriangleRTree(mesh);

    BOOST_TEST(tree->size() == 4);

    Eigen::VectorXd                        searchVector(Eigen::Vector3d(0.7, 0.5, 0.0));
    std::vector<std::pair<double, size_t>> results;

    tree->query(bgi::nearest(searchVector, 3), boost::make_function_output_iterator([&](const precice::mesh::rtree::triangle_traits::IndexType &val) {
                  results.push_back(std::make_pair(
                      boost::geometry::distance(
                          searchVector,
                          mesh->triangles()[val.second]),
                      val.second));
                }));

    std::sort(results.begin(), results.end());
    BOOST_TEST_INFO(results);
    BOOST_TEST(results.size() == 3);
    BOOST_TEST(results[0].second == tlb.getID());
    BOOST_TEST(results[1].second == tlt.getID());
    BOOST_TEST(results[2].second == trt.getID());
    BOOST_TEST(results[2].second != trb.getID());
  }
}

/// Resembles how boost geometry is used inside the PetRBF
BOOST_AUTO_TEST_CASE(QueryWithBox)
{
  PRECICE_TEST(1_rank);
  PtrMesh mesh(new precice::mesh::Mesh("MyMesh", 3, false, precice::testing::nextMeshID()));
  mesh->createVertex(Eigen::Vector3d(0, 0, 0));
  mesh->createVertex(Eigen::Vector3d(0, 0, 1));
  mesh->createVertex(Eigen::Vector3d(0, 1, 0));
  mesh->createVertex(Eigen::Vector3d(0, 1, 1));
  mesh->createVertex(Eigen::Vector3d(1, 0, 0));
  mesh->createVertex(Eigen::Vector3d(1, 0, 1));
  mesh->createVertex(Eigen::Vector3d(1, 1, 0));
  mesh->createVertex(Eigen::Vector3d(1, 1, 1));

  auto tree = rtree::getVertexRTree(mesh);
  BOOST_TEST(tree->size() == 8);

  Eigen::VectorXd searchVector(Eigen::Vector3d(0.8, 1, 0));

  {
    std::vector<size_t>             results;
    double                          radius = 0.1; // No vertices in radius
    bg::model::box<Eigen::VectorXd> search_box(
        searchVector - Eigen::VectorXd::Constant(3, radius),
        searchVector + Eigen::VectorXd::Constant(3, radius));

    tree->query(bg::index::within(search_box) and bg::index::satisfies([&](size_t const i) { return bg::distance(searchVector, mesh->vertices()[i]) <= radius; }),
                std::back_inserter(results));

    BOOST_TEST(results.empty());
  }

  {
    std::vector<size_t>             results;
    double                          radius = 0.81; // Two vertices in radius
    bg::model::box<Eigen::VectorXd> search_box(
        searchVector - Eigen::VectorXd::Constant(3, radius),
        searchVector + Eigen::VectorXd::Constant(3, radius));

    tree->query(bg::index::within(search_box) and bg::index::satisfies([&](size_t const i) { return bg::distance(searchVector, mesh->vertices()[i]) <= radius; }),
                std::back_inserter(results));

    BOOST_TEST(results.size() == 2);
    BOOST_TEST(mesh->vertices()[results[0]].getCoords() == Eigen::Vector3d(0, 1, 0));
    BOOST_TEST(mesh->vertices()[results[1]].getCoords() == Eigen::Vector3d(1, 1, 0));
  }

  {
    std::vector<size_t>             results;
    double                          radius = std::numeric_limits<double>::max();
    bg::model::box<Eigen::VectorXd> search_box(
        searchVector - Eigen::VectorXd::Constant(3, radius),
        searchVector + Eigen::VectorXd::Constant(3, radius));

    tree->query(bg::index::within(search_box) and bg::index::satisfies([&](size_t const i) { return bg::distance(searchVector, mesh->vertices()[i]) <= radius; }),
                std::back_inserter(results));

    BOOST_TEST(results.size() == 8);
  }
}

BOOST_AUTO_TEST_CASE(CacheClearing)
{
  PRECICE_TEST(1_rank);
  PtrMesh mesh(new precice::mesh::Mesh("MyMesh", 2, false, precice::testing::nextMeshID()));
  mesh->createVertex(Eigen::Vector2d(0, 0));

  // The Cache should clear whenever a mesh changes
  auto vTree1 = rtree::getVertexRTree(mesh);
  auto pTree1 = rtree::getPrimitiveRTree(mesh);
  BOOST_TEST(rtree::_cached_trees.size() == 1);
  BOOST_TEST(rtree::_primitive_trees.size() == 1);
  mesh->meshChanged(*mesh); // Emit signal, that mesh has changed
  BOOST_TEST(rtree::_cached_trees.empty());
  BOOST_TEST(rtree::_primitive_trees.empty());

  // The Cache should clear whenever we destroy the Mesh
  auto vTree2 = rtree::getVertexRTree(mesh);
  auto pTree2 = rtree::getPrimitiveRTree(mesh);
  BOOST_TEST(rtree::_cached_trees.size() == 1);
  BOOST_TEST(rtree::_primitive_trees.size() == 1);
  mesh.reset(); // Destroy mesh object, signal is emitted to clear cache
  BOOST_TEST(rtree::_cached_trees.empty());
  BOOST_TEST(rtree::_primitive_trees.empty());
}

BOOST_AUTO_TEST_CASE(PrimitveIndexComparison)
{
  PRECICE_TEST(1_rank);
  PrimitiveIndex a{Primitive::Vertex, 2lu};
  PrimitiveIndex b{Primitive::Vertex, 2lu};
  PrimitiveIndex c{Primitive::Edge, 2lu};
  PrimitiveIndex d{Primitive::Edge, 0lu};

  BOOST_TEST(a == b);
  BOOST_TEST(a != c);
  BOOST_TEST(b != c);
  BOOST_TEST(c != d);
}

BOOST_FIXTURE_TEST_CASE(IndexSinglePrimitiveType, MeshFixture)
{
  PRECICE_TEST(1_rank);
  PrimitiveRTree            rtree;
  mesh::impl::AABBGenerator gen{mesh};

  using mesh::impl::indexPrimitive;
  BOOST_TEST(rtree.empty());
  indexPrimitive(rtree, gen, mesh.vertices());
  BOOST_TEST(rtree.size() == vertex_cnt);
  indexPrimitive(rtree, gen, mesh.edges());
  BOOST_TEST(rtree.size() == vertex_cnt + edge_cnt);
  indexPrimitive(rtree, gen, mesh.triangles());
  BOOST_TEST(rtree.size() == vertex_cnt + edge_cnt + triangle_cnt);
  indexPrimitive(rtree, gen, mesh.quads());
  BOOST_TEST(rtree.size() == primitive_cnt);
}

BOOST_FIXTURE_TEST_CASE(IndexMesh, MeshFixture)
{
  PRECICE_TEST(1_rank);
  auto tree = indexMesh(mesh);
  BOOST_TEST(tree.size() == primitive_cnt);
}

BOOST_FIXTURE_TEST_CASE(CacheFunctionality, MeshFixture)
{
  PRECICE_TEST(1_rank);
  PtrMesh ptr{&mesh, [](Mesh *) {}}; // Use an empty deleter to prevent double-free

  auto vt1 = rtree::getVertexRTree(ptr);
  auto vt2 = rtree::getVertexRTree(ptr);
  BOOST_TEST(vt1 == vt2);

  auto pt1 = rtree::getPrimitiveRTree(ptr);
  auto pt2 = rtree::getPrimitiveRTree(ptr);
  BOOST_TEST(pt1 == pt2);
}

BOOST_AUTO_TEST_SUITE_END() // RTree
BOOST_AUTO_TEST_SUITE_END() // Mesh
