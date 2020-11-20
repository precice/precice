#include <Eigen/Core>
#include <algorithm>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include "io/Export.hpp"
#include "io/ExportVTK.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "query/FindClosest.hpp"
#include "query/FindClosestEdge.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

namespace precice {
namespace mesh {
class Triangle;
} // namespace mesh
} // namespace precice

using namespace precice;
using namespace precice::query;

BOOST_AUTO_TEST_SUITE(QueryTests)
BOOST_AUTO_TEST_SUITE(FindClosestTests)

BOOST_AUTO_TEST_CASE(FindClosestDistanceToVertices)
{
  PRECICE_TEST(1_rank);
  for (int dim = 2; dim <= 3; dim++) {
    mesh::Mesh mesh("RootMesh", dim, false, testing::nextMeshID());
    mesh.createVertex(Eigen::VectorXd::Zero(dim));
    Eigen::VectorXd queryCoords0 = Eigen::VectorXd::Zero(dim);
    queryCoords0(0)              = 1.0;
    FindClosest find(queryCoords0);
    find(mesh);
    query::ClosestElement closest  = find.getClosest();
    double                distance = closest.distance;
    BOOST_TEST(distance == 1.0);
    if (dim == 3) {
      mesh::Mesh mesh3D("Mesh3D", dim, false, testing::nextMeshID());
      mesh3D.createVertex(Eigen::Vector3d::Constant(1));
      FindClosest find2(Eigen::Vector3d::Constant(-1));
      find2(mesh3D);
      distance = find2.getClosest().distance;
      BOOST_TEST(distance == Eigen::Vector3d::Constant(2).norm());
    }
  }
}

BOOST_AUTO_TEST_CASE(SignOfShortestDistance)
{
  PRECICE_TEST(1_rank);
  for (int dim = 2; dim <= 3; dim++) {
    mesh::Mesh      mesh("Mesh", dim, false, testing::nextMeshID());
    mesh::Vertex &  vertex = mesh.createVertex(Eigen::VectorXd::Zero(dim));
    Eigen::VectorXd normal = Eigen::VectorXd::Zero(dim);
    normal(0)              = 1.0;
    vertex.setNormal(normal);

    // Check point that lies outside of mesh
    Eigen::VectorXd queryCoords = Eigen::VectorXd::Zero(dim);
    queryCoords(0)              = 1.0;
    FindClosest find(queryCoords);
    find(mesh);
    double distance = find.getClosest().distance;
    BOOST_TEST(distance == 1.0);

    // Check point that lies inside of mesh
    normal(0) = -1.0;
    vertex.setNormal(normal);
    find.reset();
    find(mesh);
    BOOST_TEST(find.getClosest().distance == -1.0);
  }
}

BOOST_AUTO_TEST_CASE(IndependenceOfSignOfShortestDistance)
{
  PRECICE_TEST(1_rank);
  for (int dim = 2; dim <= 3; dim++) {
    mesh::Mesh    mesh("Mesh", dim, false, testing::nextMeshID());
    mesh::Vertex &vertex = mesh.createVertex(Eigen::VectorXd::Constant(dim, 1));
    vertex.setNormal(Eigen::VectorXd::Constant(dim, 1));
    mesh::Vertex &vertex2 = mesh.createVertex(Eigen::VectorXd::Constant(dim, -2));
    vertex2.setNormal(Eigen::VectorXd::Constant(dim, 1));

    FindClosest find(Eigen::VectorXd::Zero(dim));
    find(mesh);
    double distance = find.getClosest().distance;
    BOOST_TEST(distance == -1.0 * Eigen::VectorXd::Constant(dim, 1).norm());

    // Invert normal of futher away vertex, should have no influence
    Eigen::VectorXd normal = Eigen::VectorXd::Zero(dim);
    normal(1)              = -1.0;
    vertex2.setNormal(normal);
    find.reset();
    find(mesh);
    distance = find.getClosest().distance;
    BOOST_TEST(distance == -1.0 * Eigen::VectorXd::Constant(dim, 1).norm());

    // Invert normal of closer vertex, should invert sign of distance
    vertex.setNormal(normal);
    find.reset();
    find(mesh);
    distance = find.getClosest().distance;
    BOOST_TEST(distance == Eigen::VectorXd::Constant(dim, 1).norm());
  }
}

BOOST_AUTO_TEST_CASE(FindClosestDistanceToEdges)
{
  PRECICE_TEST(1_rank);
  for (int dim = 2; dim <= 3; dim++) {
    // Create mesh consisting of two vertices and an edge
    mesh::Mesh      mesh("Mesh", dim, false, testing::nextMeshID());
    mesh::Vertex &  v1     = mesh.createVertex(Eigen::VectorXd::Constant(dim, -1));
    mesh::Vertex &  v2     = mesh.createVertex(Eigen::VectorXd::Constant(dim, 1));
    mesh::Edge &    edge   = mesh.createEdge(v1, v2);
    Eigen::VectorXd normal = Eigen::VectorXd::Constant(dim, -1);
    normal(1)              = 1.0;
    v1.setNormal(normal);
    v2.setNormal(normal);
    edge.setNormal(normal);

    // Create query points
    // Query point 0 lies outside of the mesh
    Eigen::VectorXd queryCoords0 = Eigen::VectorXd::Constant(dim, -1);
    queryCoords0(1)              = 1.0;
    // Query point 1 lies inside of the mesh
    Eigen::VectorXd queryCoords1 = Eigen::VectorXd::Constant(dim, 1);
    queryCoords1(1)              = -1.0;
    // Query point 2 lies on the query edge
    Eigen::VectorXd queryCoords2 = Eigen::VectorXd::Constant(dim, 0);
    // Query point 3 lies on the same line as the edge, but above
    Eigen::VectorXd queryCoords3 = Eigen::VectorXd::Constant(dim, 2);
    // Query point 4 lies on the same line as the edge, but below
    Eigen::VectorXd queryCoords4 = Eigen::VectorXd::Constant(dim, -2);

    // Create query objects
    FindClosest find0(queryCoords0);
    FindClosest find1(queryCoords1);
    FindClosest find2(queryCoords2);
    FindClosest find3(queryCoords3);
    FindClosest find4(queryCoords4);

    // Perform queries
    find0(mesh);
    find1(mesh);
    find2(mesh);
    find3(mesh);
    find4(mesh);

    // Evaluate query results
    Eigen::VectorXd expected = Eigen::VectorXd::Constant(dim, 1);
    BOOST_TEST(find0.getClosest().distance == expected.norm());
    BOOST_TEST(find1.getClosest().distance == -1.0 * expected.norm());
    BOOST_TEST(find2.getClosest().distance == 0.0);
    BOOST_TEST(std::abs(find3.getClosest().distance) == expected.norm());
    BOOST_TEST(std::abs(find4.getClosest().distance) == expected.norm());
  }
}

BOOST_AUTO_TEST_CASE(FindClosestDistanceToEdges3D)
{
  PRECICE_TEST(1_rank);
  int dim = 3;
  // Cremeshetry consisting of two vertices and an edge
  mesh::Mesh      mesh("Mesh", dim, false, testing::nextMeshID());
  mesh::Vertex &  v1   = mesh.createVertex(Eigen::Vector3d(-1.0, -1.0, 0.0));
  mesh::Vertex &  v2   = mesh.createVertex(Eigen::Vector3d(1.0, 1.0, 0.0));
  mesh::Edge &    edge = mesh.createEdge(v1, v2);
  Eigen::VectorXd normal(dim);
  normal << -1.0, 1.0, 0.0;
  v1.setNormal(normal);
  v2.setNormal(normal);
  edge.setNormal(normal);

  io::ExportVTK exportMesh(true);
  std::string   location = "";
  exportMesh.doExport("query-FindClosestTest-testFindClosestDistanceToEdges3D", location, mesh);

  // Create query points
  std::vector<Eigen::Vector3d> queryPoints;
  // Query point 0 lies outside of the mesh
  queryPoints.push_back(Eigen::Vector3d(-0.5, 0.0, 0.0));
  queryPoints.push_back(Eigen::Vector3d(0.0, 0.5, 0.0));
  queryPoints.push_back(Eigen::Vector3d(-0.5, 0.5, 0.0));
  queryPoints.push_back(Eigen::Vector3d(-0.5, 0.5, -0.5));
  queryPoints.push_back(Eigen::Vector3d(-0.5, 0.5, 0.5));

  // Create query objects
  std::vector<query::FindClosest *> finds;
  for (auto &queryPoint : queryPoints) {
    finds.push_back(new FindClosest(queryPoint));
  }

  // Perform queries
  for (size_t i = 0; i < finds.size(); i++) {
    BOOST_TEST((*finds.at(i))(mesh));
  }

  // Evaluate query results
  BOOST_TEST(finds.at(0)->getClosest().distance == std::sqrt(1.0 / 8.0));
  BOOST_TEST(finds.at(1)->getClosest().distance == std::sqrt(1.0 / 8.0));
  BOOST_TEST(finds.at(2)->getClosest().distance == Eigen::Vector3d(0.5, -0.5, 0.0).norm());
  BOOST_TEST(finds.at(3)->getClosest().distance == Eigen::Vector3d(0.5, -0.5, 0.5).norm());
  BOOST_TEST(finds.at(4)->getClosest().distance == Eigen::Vector3d(0.5, -0.5, -0.5).norm());

  // Clean up
  for (auto &find : finds) {
    delete find;
  }
}

BOOST_AUTO_TEST_CASE(FindClosestDistanceToTriangles)
{
  PRECICE_TEST(1_rank);
  // Create mesh to query
  mesh::Mesh    mesh("Mesh", 3, true, testing::nextMeshID());
  mesh::Vertex &v0 = mesh.createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  mesh::Vertex &v1 = mesh.createVertex(Eigen::Vector3d(1.0, 1.0, 0.0));
  mesh::Vertex &v2 = mesh.createVertex(Eigen::Vector3d(1.0, 1.0, 1.0));
  mesh::Edge &  e0 = mesh.createEdge(v0, v1);
  mesh::Edge &  e1 = mesh.createEdge(v1, v2);
  mesh::Edge &  e2 = mesh.createEdge(v2, v0);
  mesh.createTriangle(e0, e1, e2);
  mesh.computeState();

  // Prepare and issue queries
  std::vector<Eigen::Vector3d> queries;
  queries.push_back(Eigen::Vector3d(0.6, 0.6, 0.5));    //  0: on triangle middle
  queries.push_back(Eigen::Vector3d(0.0, 0.0, 0.0));    //  1: on vertex0
  queries.push_back(Eigen::Vector3d(1.0, 1.0, 0.0));    //  2: on vertex1
  queries.push_back(Eigen::Vector3d(1.0, 1.0, 1.0));    //  3: on vertex2
  queries.push_back(Eigen::Vector3d(0.5, 0.5, 0.0));    //  4: on edge0
  queries.push_back(Eigen::Vector3d(1.0, 1.0, 0.5));    //  5: on edge1
  queries.push_back(Eigen::Vector3d(0.5, 0.5, 0.5));    //  6: on edge2
  queries.push_back(Eigen::Vector3d(0.0, 1.0, 0.3));    //  7: outside triangle
  queries.push_back(Eigen::Vector3d(1.0, 0.0, 0.3));    //  8: inside triangle
  queries.push_back(Eigen::Vector3d(-0.5, -0.5, -0.5)); //  9: outside vertex0
  queries.push_back(Eigen::Vector3d(1.5, 1.5, -0.5));   // 10: outside vertex1
  queries.push_back(Eigen::Vector3d(1.5, 1.5, 1.5));    // 11: outside vertex2
  std::vector<std::shared_ptr<FindClosest>> findVisitors;
  for (size_t i = 0; i < queries.size(); i++) {
    std::shared_ptr<FindClosest> find(new FindClosest(queries.at(i)));
    findVisitors.push_back(find);
    (*findVisitors.at(i))(mesh);
  }

  // Validate results
  for (size_t i = 0; i < 7; i++) {
    BOOST_TEST(findVisitors.at(i)->getClosest().distance == 0.0);
  }
  Eigen::Vector3d expect(0.5, -0.5, 0.0);
  BOOST_TEST(testing::equals(findVisitors.at(7)->getClosest().vectorToElement, expect));
  BOOST_TEST(findVisitors.at(7)->getClosest().distance == expect.norm());
  expect << -0.5, 0.5, 0.0;
  BOOST_TEST(testing::equals(findVisitors.at(8)->getClosest().vectorToElement, expect));
  BOOST_TEST(findVisitors.at(8)->getClosest().distance == -1.0 * expect.norm());
  expect << 0.5, 0.5, 0.5;
  BOOST_TEST(testing::equals(findVisitors.at(9)->getClosest().vectorToElement, expect));
  BOOST_TEST(std::abs(findVisitors.at(9)->getClosest().distance) == expect.norm());
  expect << -0.5, -0.5, 0.5;
  BOOST_TEST(testing::equals(findVisitors.at(10)->getClosest().vectorToElement, expect));
  BOOST_TEST(std::abs(findVisitors.at(10)->getClosest().distance) == expect.norm());
  expect << -0.5, -0.5, -0.5;
  BOOST_TEST(testing::equals(findVisitors.at(11)->getClosest().vectorToElement, expect));
  BOOST_TEST(std::abs(findVisitors.at(11)->getClosest().distance) == expect.norm());
}

BOOST_AUTO_TEST_CASE(FindClosestDistanceToTrianglesAndVertices)
{
  PRECICE_TEST(1_rank);
  int           dim = 2;
  mesh::Mesh    mesh("Mesh", dim, false, testing::nextMeshID());
  mesh::Vertex &vertex1 = mesh.createVertex(Eigen::Vector2d(0.0, 0.0));
  vertex1.setNormal(Eigen::Vector2d(-0.5, 0.5));

  mesh::Vertex &vertex2 = mesh.createVertex(Eigen::Vector2d(1.0, 0.0));
  vertex2.setNormal(Eigen::Vector2d(0.5, 0.5));

  mesh::Edge &edge = mesh.createEdge(vertex1, vertex2);
  edge.setNormal(Eigen::Vector2d(0.0, 1.0));

  query::FindClosest find(Eigen::Vector2d(0.0, 0.0));
  find(mesh);
  double distance = find.getClosest().distance;
  BOOST_TEST(distance == 0.0);

  query::FindClosest find2(Eigen::Vector2d(0.5, 0.0));
  find2(mesh);
  distance = find2.getClosest().distance;
  BOOST_TEST(distance == 0.0);

  query::FindClosest find3(Eigen::Vector2d(0.5, 0.1));
  find3(mesh);
  distance = find3.getClosest().distance;
  BOOST_TEST(distance == 0.1);

  query::FindClosest find4(Eigen::Vector2d(0.0, 1.5));
  find4(mesh);
  distance = find4.getClosest().distance;
  BOOST_TEST(distance == 1.5);

  query::FindClosest find5(Eigen::Vector2d(0.5, -1.0));
  find5(mesh);
  distance = find5.getClosest().distance;
  BOOST_TEST(distance == -1.0);
}

BOOST_AUTO_TEST_CASE(WeigthsOfVertices)
{
  PRECICE_TEST(1_rank);
  int dim = 2;

  // Create mesh
  mesh::Mesh    mesh("Mesh", dim, true, testing::nextMeshID());
  mesh::Vertex &vertex1 = mesh.createVertex(Eigen::Vector2d(0.0, 0.0));
  mesh::Vertex &vertex2 = mesh.createVertex(Eigen::Vector2d(1.0, 0.0));
  mesh.createEdge(vertex1, vertex2);
  mesh.computeState();

  // Query elements
  query::FindClosest findClosest(Eigen::Vector2d(0.3, 1.0));
  findClosest(mesh);
  const query::ClosestElement &closest = findClosest.getClosest();

  // Validate results
  BOOST_TEST(closest.interpolationElements.size() == 2);
  auto pointerVertex1 = closest.interpolationElements.at(0).element;
  auto pointerVertex2 = closest.interpolationElements.at(1).element;
  BOOST_TEST(testing::equals(pointerVertex1->getCoords(), Eigen::Vector2d(0.0, 0.0)));
  BOOST_TEST(testing::equals(pointerVertex2->getCoords(), Eigen::Vector2d(1.0, 0.0)));
  BOOST_TEST(closest.interpolationElements.at(0).weight == 0.7);
  BOOST_TEST(closest.interpolationElements.at(1).weight == 0.3);
}

struct MeshFixture {
  int             dimension = 3;
  double          z         = 0.0;
  mesh::Mesh      mesh;
  mesh::Vertex *  v1, *v2, *v3, *vinside, *voutside;
  mesh::Edge *    e12, *e23, *e31;
  mesh::Triangle *t;
  MeshFixture()
      : mesh("Mesh", dimension, true, testing::nextMeshID())
  {
    v1       = &mesh.createVertex(Eigen::Vector3d(0.0, 0.0, z));
    v2       = &mesh.createVertex(Eigen::Vector3d(1.0, 0.0, z));
    v3       = &mesh.createVertex(Eigen::Vector3d(0.5, 0.5, z));
    vinside  = &mesh.createVertex(Eigen::Vector3d(0.1, 0.1, z));
    voutside = &mesh.createVertex(Eigen::Vector3d(1.0, 1.0, z));
    e12      = &mesh.createEdge(*v1, *v2);
    e23      = &mesh.createEdge(*v2, *v3);
    e31      = &mesh.createEdge(*v3, *v1);
    t        = &mesh.createTriangle(*e12, *e23, *e31);
    mesh.computeState();
  }
  ~MeshFixture() {}
};

BOOST_FIXTURE_TEST_SUITE(InterpolationElements, MeshFixture)

BOOST_AUTO_TEST_CASE(Vertex)
{
  PRECICE_TEST(1_rank);
  auto elems = query::generateInterpolationElements(*vinside, *v1);
  BOOST_TEST(elems.size() == 1);
  BOOST_TEST(elems.front().element == v1);
  BOOST_TEST(elems.front().weight == 1.0);
}

BOOST_AUTO_TEST_CASE(Edge)
{
  PRECICE_TEST(1_rank);
  auto elems = query::generateInterpolationElements(*vinside, *e12);
  BOOST_TEST(elems.size() == 2);
  BOOST_TEST(elems.front().element == v1);
  BOOST_TEST(elems.front().weight == 0.9);
  BOOST_TEST(elems.back().element == v2);
  BOOST_TEST(elems.back().weight == 0.1);
}

BOOST_AUTO_TEST_CASE(TriangleInside)
{
  PRECICE_TEST(1_rank);
  auto elems = query::generateInterpolationElements(*vinside, *t);
  BOOST_TEST(elems.size() == 3);
  std::map<mesh::Vertex const *, double> weights;
  for (const auto &elem : elems) {
    weights[elem.element] = elem.weight;
  }
  BOOST_TEST(weights.at(v1) == 0.8);
  BOOST_TEST(weights.at(v2) == 0.0);
  BOOST_TEST(weights.at(v3) == 0.2);
}

BOOST_AUTO_TEST_CASE(TriangleOutside)
{
  PRECICE_TEST(1_rank);
  auto elems = query::generateInterpolationElements(*voutside, *t);
  BOOST_TEST(elems.size() == 3);
  std::map<mesh::Vertex const *, double> weights;
  for (const auto &elem : elems) {
    weights[elem.element] = elem.weight;
  }
  // Extrapolating
  BOOST_TEST(weights.at(v1) == -1.0);
  BOOST_TEST(weights.at(v2) == 0.0);
  BOOST_TEST(weights.at(v3) == 2.0);
}

BOOST_AUTO_TEST_SUITE_END() // InterpolationElements

BOOST_AUTO_TEST_SUITE_END() // FindClosestTests
BOOST_AUTO_TEST_SUITE_END() // QueryTests
