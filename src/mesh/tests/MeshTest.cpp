#include <Eigen/Core>
#include <algorithm>
#include <deque>
#include <iosfwd>
#include <memory>
#include <string>
#include <vector>
#include "logging/Logger.hpp"
#include "mesh/BoundingBox.hpp"
#include "mesh/Data.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Quad.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Vertex.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

using namespace precice;
using namespace precice::mesh;
using Eigen::Vector2d;
using Eigen::Vector3d;
using precice::testing::equals;

BOOST_AUTO_TEST_SUITE(MeshTests)

BOOST_AUTO_TEST_SUITE(MeshTests)

BOOST_AUTO_TEST_CASE(ComputeState_2D)
{
  PRECICE_TEST(1_rank);
  mesh::Mesh mesh("MyMesh", 2, true, testing::nextMeshID());
  // Create mesh
  Vertex &v1 = mesh.createVertex(Vector2d(0.0, 0.0));
  Vertex &v2 = mesh.createVertex(Vector2d(1.0, 0.0));
  Vertex &v3 = mesh.createVertex(Vector2d(1.0, 1.0));
  //
  //
  // *****
  Edge &e1 = mesh.createEdge(v1, v2);
  //     *
  //     *  <---
  // *****
  Edge &e2 = mesh.createEdge(v2, v3);
  mesh.computeState();

  // Perform test validations
  BOOST_TEST(equals(e1.getCenter(), Vector2d(0.5, 0.0)));
  BOOST_TEST(equals(e2.getCenter(), Vector2d(1.0, 0.5)));
  BOOST_TEST(e1.getEnclosingRadius() == 0.5);
  BOOST_TEST(e2.getEnclosingRadius() == 0.5);
  BOOST_TEST(equals(e1.getNormal(), Vector2d(0.0, 1.0)));
  BOOST_TEST(equals(e2.getNormal(), Vector2d(-1.0, 0.0)));
  BOOST_TEST(equals(v1.getNormal(), Vector2d(0.0, 1.0)));
  BOOST_TEST(equals(v2.getNormal(), Vector2d(-std::sqrt(0.5), std::sqrt(0.5))));
  BOOST_TEST(equals(v3.getNormal(), Vector2d(-1.0, 0.0)));
}

BOOST_AUTO_TEST_CASE(ComputeState_3D_Triangle)
{
  PRECICE_TEST(1_rank);
  precice::mesh::Mesh mesh("MyMesh", 3, true, testing::nextMeshID());
  // Create mesh
  Vertex &v1 = mesh.createVertex(Vector3d(0.0, 0.0, 0.0));
  Vertex &v2 = mesh.createVertex(Vector3d(1.0, 0.0, 1.0));
  Vertex &v3 = mesh.createVertex(Vector3d(1.0, 1.0, 1.0));
  Vertex &v4 = mesh.createVertex(Vector3d(2.0, 0.0, 2.0));
  Edge &  e1 = mesh.createEdge(v1, v2);
  Edge &  e2 = mesh.createEdge(v2, v3);
  Edge &  e3 = mesh.createEdge(v3, v1);
  Edge &  e4 = mesh.createEdge(v2, v4);
  Edge &  e5 = mesh.createEdge(v4, v3);

  //       *
  //     * *
  //   *   *
  // *******
  Triangle &t1 = mesh.createTriangle(e1, e2, e3);
  //       *
  //     * * *     <---
  //   *   *   *
  // *************
  Triangle &t2 = mesh.createTriangle(e4, e5, e2);
  mesh.computeState();

  // Perform test validations
  BOOST_TEST(equals(e1.getCenter(), Vector3d(0.5, 0.0, 0.5)));
  BOOST_TEST(equals(e2.getCenter(), Vector3d(1.0, 0.5, 1.0)));
  BOOST_TEST(equals(e3.getCenter(), Vector3d(0.5, 0.5, 0.5)));
  BOOST_TEST(equals(e4.getCenter(), Vector3d(1.5, 0.0, 1.5)));
  BOOST_TEST(equals(e5.getCenter(), Vector3d(1.5, 0.5, 1.5)));
  BOOST_TEST(e1.getEnclosingRadius() == std::sqrt(2.0) * 0.5);
  BOOST_TEST(e2.getEnclosingRadius() == 0.5);
  BOOST_TEST(e3.getEnclosingRadius() == std::sqrt(3.0) * 0.5);
  BOOST_TEST(e4.getEnclosingRadius() == std::sqrt(2.0) * 0.5);
  BOOST_TEST(e5.getEnclosingRadius() == std::sqrt(3.0) * 0.5);

  BOOST_TEST(equals(t1.getCenter(), Vector3d(2.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0)));
  BOOST_TEST(equals(t2.getCenter(), Vector3d(4.0 / 3.0, 1.0 / 3.0, 4.0 / 3.0)));
  BOOST_TEST(t1.getEnclosingRadius() == 1.0);
  BOOST_TEST(t2.getEnclosingRadius() == 1.0);
  Vector3d normal(1.0, 0.0, -1.0);
  normal = normal.normalized();
  BOOST_TEST(normal.norm() == 1.0);
  BOOST_TEST(equals(t1.getNormal(), normal));
  BOOST_TEST(equals(t2.getNormal(), normal));

  BOOST_TEST(equals(e1.getNormal(), normal));
  BOOST_TEST(equals(e2.getNormal(), normal));
  BOOST_TEST(equals(e3.getNormal(), normal));
  BOOST_TEST(equals(e4.getNormal(), normal));
  BOOST_TEST(equals(e5.getNormal(), normal));
  BOOST_TEST(equals(v1.getNormal(), normal));
  BOOST_TEST(equals(v2.getNormal(), normal));
  BOOST_TEST(equals(v3.getNormal(), normal));
  BOOST_TEST(equals(v4.getNormal(), normal));
}

BOOST_AUTO_TEST_CASE(ComputeState_3D_Quad)
{
  PRECICE_TEST(1_rank);
  mesh::Mesh mesh("MyMesh", 3, true, testing::nextMeshID());
  // Create mesh (Two rectangles with a common edge at z-axis. One extends in
  // x-y-plane (2 long), the other in y-z-plane (1 long))
  Vertex &v0 = mesh.createVertex(Vector3d(0.0, 0.0, 0.0));
  Vertex &v1 = mesh.createVertex(Vector3d(2.0, 0.0, 0.0));
  Vertex &v2 = mesh.createVertex(Vector3d(2.0, 1.0, 0.0));
  Vertex &v3 = mesh.createVertex(Vector3d(0.0, 1.0, 0.0));
  Vertex &v4 = mesh.createVertex(Vector3d(0.0, 0.0, 1.0));
  Vertex &v5 = mesh.createVertex(Vector3d(0.0, 1.0, 1.0));
  Edge &  e0 = mesh.createEdge(v0, v1);
  Edge &  e1 = mesh.createEdge(v1, v2);
  Edge &  e2 = mesh.createEdge(v2, v3);
  Edge &  e3 = mesh.createEdge(v3, v0);
  Edge &  e4 = mesh.createEdge(v0, v4);
  Edge &  e5 = mesh.createEdge(v4, v5);
  Edge &  e6 = mesh.createEdge(v5, v3);

  Quad &q0 = mesh.createQuad(e0, e1, e2, e3); // in x-y-plane
  Quad &q1 = mesh.createQuad(e4, e5, e6, e3); // in z-y-plane
  mesh.computeState();

  // Perform test validations
  BOOST_TEST(equals(e0.getCenter(), Vector3d(1.0, 0.0, 0.0)));
  BOOST_TEST(equals(e1.getCenter(), Vector3d(2.0, 0.5, 0.0)));
  BOOST_TEST(equals(e2.getCenter(), Vector3d(1.0, 1.0, 0.0)));
  BOOST_TEST(equals(e3.getCenter(), Vector3d(0.0, 0.5, 0.0)));
  BOOST_TEST(equals(e4.getCenter(), Vector3d(0.0, 0.0, 0.5)));
  BOOST_TEST(equals(e5.getCenter(), Vector3d(0.0, 0.5, 1.0)));
  BOOST_TEST(equals(e6.getCenter(), Vector3d(0.0, 1.0, 0.5)));

  BOOST_TEST(e0.getEnclosingRadius() == 1.0);
  BOOST_TEST(e1.getEnclosingRadius() == 0.5);
  BOOST_TEST(e2.getEnclosingRadius() == 1.0);
  BOOST_TEST(e3.getEnclosingRadius() == 0.5);
  BOOST_TEST(e4.getEnclosingRadius() == 0.5);
  BOOST_TEST(e5.getEnclosingRadius() == 0.5);
  BOOST_TEST(e6.getEnclosingRadius() == 0.5);

  BOOST_TEST(equals(q0.getCenter(), Vector3d(1.0, 0.5, 0.0)));
  BOOST_TEST(equals(q1.getCenter(), Vector3d(0.0, 0.5, 0.5)));

  BOOST_TEST(q0.getEnclosingRadius() == sqrt(1.25));
  BOOST_TEST(q1.getEnclosingRadius() == sqrt(0.5));

  Vector3d normal0(0.0, 0.0, -1.0);
  BOOST_TEST(equals(q0.getNormal(), normal0));
  Vector3d normal1(1.0, 0.0, 0.0);
  BOOST_TEST(equals(q1.getNormal(), normal1));

  BOOST_TEST(equals(e0.getNormal(), normal0));
  BOOST_TEST(equals(e1.getNormal(), normal0));
  BOOST_TEST(equals(e2.getNormal(), normal0));
  Vector3d normal0and1(2.0 * normal0 + normal1);
  normal0and1 = normal0and1.normalized();
  BOOST_TEST(equals(e3.getNormal(), normal0and1));
  BOOST_TEST(equals(e4.getNormal(), normal1));
  BOOST_TEST(equals(e5.getNormal(), normal1));
  BOOST_TEST(equals(e6.getNormal(), normal1));

  BOOST_TEST(equals(v0.getNormal(), normal0and1));
  BOOST_TEST(equals(v1.getNormal(), normal0));
  BOOST_TEST(equals(v2.getNormal(), normal0));
  BOOST_TEST(equals(v3.getNormal(), normal0and1));
  BOOST_TEST(equals(v4.getNormal(), normal1));
  BOOST_TEST(equals(v5.getNormal(), normal1));
}

BOOST_AUTO_TEST_CASE(BoundingBoxCOG_2D)
{
  PRECICE_TEST(1_rank);
  Eigen::Vector2d coords0(2, 0);
  Eigen::Vector2d coords1(-1, 4);
  Eigen::Vector2d coords2(0, 1);

  mesh::Mesh mesh("2D Testmesh", 2, false, testing::nextMeshID());
  mesh.createVertex(coords0);
  mesh.createVertex(coords1);
  mesh.createVertex(coords2);

  mesh.computeBoundingBox();

  mesh::BoundingBox bBox = mesh.getBoundingBox();
  auto              cog  = bBox.center();

  mesh::BoundingBox referenceBox({-1.0, 2.0,
                                  0.0, 4.0});

  std::vector<double> referenceCOG = {0.5, 2.0};

  BOOST_TEST(bBox.getDimension() == 2);
  BOOST_TEST(cog.size() == 2);
  BOOST_TEST(referenceBox == bBox);

  for (decltype(cog.size()) d = 0; d < cog.size(); d++) {
    BOOST_TEST(referenceCOG[d] == cog[d]);
  }
}

BOOST_AUTO_TEST_CASE(BoundingBoxCOG_3D)
{
  PRECICE_TEST(1_rank);
  Eigen::Vector3d coords0(2, 0, -3);
  Eigen::Vector3d coords1(-1, 4, 8);
  Eigen::Vector3d coords2(0, 1, -2);
  Eigen::Vector3d coords3(3.5, 2, -2);

  mesh::Mesh mesh("3D Testmesh", 3, false, testing::nextMeshID());
  mesh.createVertex(coords0);
  mesh.createVertex(coords1);
  mesh.createVertex(coords2);
  mesh.createVertex(coords3);

  mesh.computeBoundingBox();

  mesh::BoundingBox bBox = mesh.getBoundingBox();
  auto              cog  = bBox.center();

  mesh::BoundingBox referenceBox({-1.0, 3.5,
                                  0.0, 4.0,
                                  -3.0, 8.0});

  std::vector<double> referenceCOG = {1.25, 2.0, 2.5};

  BOOST_TEST(bBox.getDimension() == 3);
  BOOST_TEST(cog.size() == 3);
  BOOST_TEST(referenceBox == bBox);

  for (decltype(cog.size()) d = 0; d < cog.size(); d++) {
    BOOST_TEST(referenceCOG[d] == cog[d]);
  }
}

BOOST_AUTO_TEST_CASE(Demonstration)
{
  PRECICE_TEST(1_rank);
  for (int dim = 2; dim <= 3; dim++) {
    // Create mesh object
    std::string         meshName("MyMesh");
    bool                flipNormals = false; // The normals of triangles, edges, vertices
    precice::mesh::Mesh mesh(meshName, dim, flipNormals, testing::nextMeshID());

    // Validate mesh object state
    BOOST_TEST(mesh.getName() == meshName);
    BOOST_TEST(mesh.isFlipNormals() == flipNormals);

    // Create mesh vertices
    Eigen::VectorXd coords0(dim);
    Eigen::VectorXd coords1(dim);
    Eigen::VectorXd coords2(dim);
    if (dim == 2) {
      coords0 << 0.0, 0.0;
      coords1 << 1.0, 0.0;
      coords2 << 0.0, 1.0;
    } else {
      coords0 << 0.0, 0.0, 0.0;
      coords1 << 1.0, 0.0, 0.0;
      coords2 << 0.0, 0.0, 1.0;
    }
    Vertex &v0 = mesh.createVertex(coords0);
    Vertex &v1 = mesh.createVertex(coords1);
    Vertex &v2 = mesh.createVertex(coords2);

    // Validate mesh vertices state
    // This is the preferred way to iterate over elements in a mesh, it hides
    // the details of the vertex container in class Mesh.
    size_t index = 0;
    for (Vertex &vertex : mesh.vertices()) {
      if (index == 0) {
        BOOST_TEST(vertex.getID() == v0.getID());
      } else if (index == 1) {
        BOOST_TEST(vertex.getID() == v1.getID());
      } else if (index == 2) {
        BOOST_TEST(vertex.getID() == v2.getID());
      } else {
        BOOST_TEST(false);
      }
      index++;
    }

    // Create mesh edges
    Edge &e0 = mesh.createEdge(v0, v1);
    Edge &e1 = mesh.createEdge(v1, v2);
    Edge &e2 = mesh.createEdge(v2, v0);

    // Validate mesh edges state
    index = 0;
    for (Edge &edge : mesh.edges()) {
      if (index == 0) {
        BOOST_TEST(edge.getID() == e0.getID());
      } else if (index == 1) {
        BOOST_TEST(edge.getID() == e1.getID());
      } else if (index == 2) {
        BOOST_TEST(edge.getID() == e2.getID());
      } else {
        BOOST_TEST(false);
      }
      index++;
    }

    Triangle *t = nullptr;
    if (dim == 3) {
      // Create triangle
      t = &mesh.createTriangle(e0, e1, e2);

      // Validate mesh triangle
      BOOST_TEST((*mesh.triangles().begin()).getID() == t->getID());
    }

    // Create vertex data
    std::string dataName("MyData");
    int         dataDimensions = dim;
    // Add a data set to the mesh. Every data value is associated to a vertex in
    // the mesh via the vertex ID. The data values are created when
    // Mesh::computeState() is called by the mesh holding the data.
    PtrData data = mesh.createData(dataName, dataDimensions);

    // Validate data state
    BOOST_TEST(data->getName() == dataName);
    BOOST_TEST(data->getDimensions() == dataDimensions);

    // Validate state of mesh with data
    BOOST_TEST(mesh.data().size() == 1);
    BOOST_TEST(mesh.data()[0]->getName() == dataName);

    // Compute the state of the mesh elements (vertices, edges, triangles)
    mesh.computeState();

    // Allocate memory for the data values of set data. Before data value access
    // leads to assertions.
    mesh.allocateDataValues();

    // Access data values
    Eigen::VectorXd &dataValues = data->values();
    BOOST_TEST(dataValues.size() == 3 * dim);
    BOOST_TEST(v0.getID() == 0);
    Eigen::VectorXd value = Eigen::VectorXd::Zero(dim);
    for (int i = 0; i < dim; i++) {
      value(i) = dataValues(v0.getID() * dim + i);
    }
  }
}

BOOST_AUTO_TEST_CASE(MeshEquality)
{
  PRECICE_TEST(1_rank);
  int   dim = 3;
  Mesh  mesh1("Mesh1", dim, false, testing::nextMeshID());
  Mesh  mesh1flipped("Mesh1flipped", dim, true, testing::nextMeshID());
  Mesh  mesh2("Mesh2", dim, false, testing::nextMeshID());
  Mesh *meshes[3] = {&mesh1, &mesh1flipped, &mesh2};
  for (auto ptr : meshes) {
    auto &          mesh = *ptr;
    Eigen::VectorXd coords0(dim);
    Eigen::VectorXd coords1(dim);
    Eigen::VectorXd coords2(dim);
    Eigen::VectorXd coords3(dim);
    coords0 << 0.0, 0.0, 0.0;
    coords1 << 1.0, 0.0, 0.0;
    coords2 << 0.0, 0.0, 1.0;
    coords3 << 1.0, 0.0, 1.0;
    Vertex &v0 = mesh.createVertex(coords0);
    Vertex &v1 = mesh.createVertex(coords1);
    Vertex &v2 = mesh.createVertex(coords2);
    Vertex &v3 = mesh.createVertex(coords3);
    Edge &  e0 = mesh.createEdge(v0, v1); // LINESTRING (0 0 0, 1 0 0)
    Edge &  e1 = mesh.createEdge(v1, v2); // LINESTRING (1 0 0, 0 0 1)
    Edge &  e2 = mesh.createEdge(v2, v0); // LINESTRING (0 0 1, 0 0 0)
    Edge &  e3 = mesh.createEdge(v1, v3); // LINESTRING (1 0 0, 1 0 1)
    Edge &  e4 = mesh.createEdge(v3, v2); // LINESTRING (1 0 1, 0 0 1)
    mesh.createTriangle(e0, e1, e2);
    mesh.createQuad(e0, e3, e4, e2);
    mesh.computeState();
  }
  BOOST_TEST(mesh1 != mesh1flipped);
  BOOST_TEST(mesh1 == mesh2);
}

BOOST_AUTO_TEST_CASE(MeshWKTPrint)
{
  PRECICE_TEST(1_rank);
  Mesh    mesh("WKTMesh", 3, false, testing::nextMeshID());
  Vertex &v0 = mesh.createVertex(Eigen::Vector3d(0., 0., 0.));
  Vertex &v1 = mesh.createVertex(Eigen::Vector3d(1., 0., 0.));
  Vertex &v2 = mesh.createVertex(Eigen::Vector3d(0., 0., 1.));
  Vertex &v3 = mesh.createVertex(Eigen::Vector3d(1., 0., 1.));
  Edge &  e0 = mesh.createEdge(v0, v1); // LINESTRING (0 0 0, 1 0 0)
  Edge &  e1 = mesh.createEdge(v1, v2); // LINESTRING (1 0 0, 0 0 1)
  Edge &  e2 = mesh.createEdge(v2, v0); // LINESTRING (0 0 1, 0 0 0)
  Edge &  e3 = mesh.createEdge(v1, v3); // LINESTRING (1 0 0, 1 0 1)
  Edge &  e4 = mesh.createEdge(v3, v2); // LINESTRING (1 0 1, 0 0 1)
  mesh.createTriangle(e0, e1, e2);
  mesh.createQuad(e0, e3, e4, e2);
  mesh.computeState();
  std::stringstream sstream;
  sstream << mesh;
  std::string reference(
      "Mesh \"WKTMesh\", dimensionality = 3:\n"
      "GEOMETRYCOLLECTION(\n"
      "POINT (0 0 0), POINT (1 0 0), POINT (0 0 1), POINT (1 0 1),\n"
      "LINESTRING (0 0 0, 1 0 0), LINESTRING (1 0 0, 0 0 1), LINESTRING (0 0 1, 0 0 0), LINESTRING (1 0 0, 1 0 1), LINESTRING (1 0 1, 0 0 1),\n"
      "POLYGON ((0 0 0, 1 0 0, 0 0 1, 0 0 0)),\n"
      "POLYGON ((0 0 0, 1 0 0, 1 0 1, 0 0 1, 0 0 0))\n"
      ")");
  BOOST_TEST(reference == sstream.str());
}

BOOST_AUTO_TEST_CASE(CreateUniqueEdge)
{
  PRECICE_TEST(1_rank);
  int             dim = 3;
  Mesh            mesh1("Mesh1", dim, false, testing::nextMeshID());
  auto &          mesh = mesh1;
  Eigen::VectorXd coords0(dim);
  Eigen::VectorXd coords1(dim);
  Eigen::VectorXd coords2(dim);
  coords0 << 0.0, 0.0, 0.0;
  coords1 << 1.0, 0.0, 0.0;
  coords2 << 0.0, 0.0, 1.0;
  Vertex &v0 = mesh.createVertex(coords0);
  Vertex &v1 = mesh.createVertex(coords1);
  Vertex &v2 = mesh.createVertex(coords2);

  Edge &e01a = mesh.createEdge(v0, v1); // LINESTRING (0 0 0, 1 0 0)
  mesh.createEdge(v0, v1);              // LINESTRING (0 0 0, 1 0 0)
  BOOST_TEST(mesh.edges().size() == 2);

  Edge &e01c = mesh.createUniqueEdge(v0, v1); // LINESTRING (0 0 0, 1 0 0)
  BOOST_TEST(mesh.edges().size() == 2);
  BOOST_TEST(e01a == e01c);

  mesh.createUniqueEdge(v1, v2); // LINESTRING (0 0 0, 1 0 0)
  BOOST_TEST(mesh.edges().size() == 3);
  mesh.createUniqueEdge(v1, v2); // LINESTRING (0 0 0, 1 0 0)
  BOOST_TEST(mesh.edges().size() == 3);
}

BOOST_AUTO_TEST_CASE(ComputeStateOfNotFullyConnectedMesh)
{
  Mesh mesh("Mesh1", 3, false, testing::nextMeshID());
  mesh.createData("Data", 1);
  const double ctz = 0.0000013;

  Eigen::Vector3d coords0;
  Eigen::Vector3d coords1;
  Eigen::Vector3d coords2;
  Eigen::Vector3d coords3;
  Eigen::Vector3d coords4;
  Eigen::Vector3d coords5;
  coords0 << ctz, ctz, ctz;
  coords1 << 1.0, ctz, ctz;
  coords2 << ctz, 1.0, ctz;
  coords3 << ctz, -1.0, ctz;
  coords4 << ctz, ctz, ctz - 1.;    // edge only
  coords5 << ctz, ctz - 1, ctz - 1; // disconnected
  Vertex &v0 = mesh.createVertex(coords0);
  Vertex &v1 = mesh.createVertex(coords1);
  Vertex &v2 = mesh.createVertex(coords2);
  Vertex &v3 = mesh.createVertex(coords3);
  Vertex &v4 = mesh.createVertex(coords4);
  mesh.createVertex(coords5);
  BOOST_TEST(mesh.vertices().size() == 6);

  Edge &e0 = mesh.createEdge(v0, v1);
  Edge &e1 = mesh.createEdge(v1, v2);
  Edge &e2 = mesh.createEdge(v2, v0);
  Edge &e3 = mesh.createEdge(v1, v3);
  Edge &e4 = mesh.createEdge(v3, v0);
  mesh.createEdge(v0, v4); // edge only
  BOOST_TEST(mesh.edges().size() == 6);

  mesh.createTriangle(e0, e1, e2);
  mesh.createTriangle(e0, e3, e4);
  BOOST_TEST(mesh.triangles().size() == 2);

  mesh.allocateDataValues();
  BOOST_TEST(mesh.data().size() == 1);
  mesh.computeState();

  for (const auto &vertex : mesh.vertices()) {
    BOOST_TEST(vertex.getNormal().allFinite());
  }
}

BOOST_AUTO_TEST_CASE(ResizeDataGrow)
{
  PRECICE_TEST(1_rank);
  precice::mesh::Mesh mesh("MyMesh", 3, true, testing::nextMeshID());
  const auto &        values = mesh.createData("Data", 1)->values();

  // Create mesh
  mesh.createVertex(Vector3d(0.0, 0.0, 0.0));
  mesh.createVertex(Vector3d(1.0, 0.0, 1.0));

  BOOST_TEST(mesh.vertices().size() == 2);
  mesh.allocateDataValues();
  BOOST_TEST(values.size() == 2);

  mesh.createVertex(Vector3d(1.0, 1.0, 1.0));
  mesh.createVertex(Vector3d(2.0, 0.0, 2.0));
  mesh.createVertex(Vector3d(2.0, 0.0, 2.1));

  BOOST_TEST(mesh.vertices().size() == 5);
  mesh.allocateDataValues();
  BOOST_TEST(values.size() == 5);
}

BOOST_AUTO_TEST_CASE(ResizeDataShrink)
{
  PRECICE_TEST(1_rank);
  precice::mesh::Mesh mesh("MyMesh", 3, true, testing::nextMeshID());
  const auto &        values = mesh.createData("Data", 1)->values();

  // Create mesh
  mesh.createVertex(Vector3d(0.0, 0.0, 0.0));
  mesh.createVertex(Vector3d(1.0, 0.0, 1.0));
  mesh.createVertex(Vector3d(1.0, 1.0, 1.0));
  mesh.createVertex(Vector3d(2.0, 2.0, 2.0));

  BOOST_TEST(mesh.vertices().size() == 4);
  mesh.allocateDataValues();
  BOOST_TEST(values.size() == 4);

  mesh.clear();
  mesh.createVertex(Vector3d(0.0, 0.0, 0.0));
  mesh.createVertex(Vector3d(1.0, 0.0, 1.0));

  BOOST_TEST(mesh.vertices().size() == 2);
  mesh.allocateDataValues();
  BOOST_TEST(values.size() == 2);
}

BOOST_AUTO_TEST_SUITE_END() // Mesh
BOOST_AUTO_TEST_SUITE_END() // Mesh
