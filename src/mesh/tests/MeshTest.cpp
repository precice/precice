#include <Eigen/Core>
#include <Eigen/src/Core/Matrix.h>
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
#include "mesh/SharedPointer.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Utils.hpp"
#include "mesh/Vertex.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "utils/algorithm.hpp"

using namespace precice;
using namespace precice::mesh;
using Eigen::Vector2d;
using Eigen::Vector3d;
using precice::testing::equals;

BOOST_AUTO_TEST_SUITE(MeshTests)

BOOST_AUTO_TEST_SUITE(MeshTests)

BOOST_AUTO_TEST_CASE(BoundingBoxCOG_2D)
{
  PRECICE_TEST(1_rank);
  Eigen::Vector2d coords0(2, 0);
  Eigen::Vector2d coords1(-1, 4);
  Eigen::Vector2d coords2(0, 1);

  mesh::Mesh mesh("2D Testmesh", 2, testing::nextMeshID());
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
    BOOST_TEST(referenceCOG.at(d) == cog(d));
  }
}

BOOST_AUTO_TEST_CASE(BoundingBoxCOG_3D)
{
  PRECICE_TEST(1_rank);
  Eigen::Vector3d coords0(2, 0, -3);
  Eigen::Vector3d coords1(-1, 4, 8);
  Eigen::Vector3d coords2(0, 1, -2);
  Eigen::Vector3d coords3(3.5, 2, -2);

  mesh::Mesh mesh("3D Testmesh", 3, testing::nextMeshID());
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
    BOOST_TEST(referenceCOG.at(d) == cog(d));
  }
}

BOOST_AUTO_TEST_CASE(Demonstration)
{
  PRECICE_TEST(1_rank);
  for (int dim = 2; dim <= 3; dim++) {
    // Create mesh object
    std::string         meshName("MyMesh");
    precice::mesh::Mesh mesh(meshName, dim, testing::nextMeshID());

    // Validate mesh object state
    BOOST_TEST(mesh.getName() == meshName);

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
    // the mesh via the vertex ID.
    PtrData data = mesh.createData(dataName, dataDimensions, 0_dataID);

    // Validate data state
    BOOST_TEST(data->getName() == dataName);
    BOOST_TEST(data->getDimensions() == dataDimensions);

    // Validate state of mesh with data
    BOOST_TEST(mesh.data().size() == 1);
    BOOST_TEST(mesh.data(0)->getName() == dataName);

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
  int                   dim = 3;
  Mesh                  mesh1("Mesh1", dim, testing::nextMeshID());
  Mesh                  mesh2("Mesh2", dim, testing::nextMeshID());
  std::array<Mesh *, 2> meshes = {&mesh1, &mesh2};
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
    mesh.createEdge(v1, v3);              // LINESTRING (1 0 0, 1 0 1)
    mesh.createEdge(v3, v2);              // LINESTRING (1 0 1, 0 0 1)
    mesh.createTriangle(e0, e1, e2);
  }
  BOOST_TEST(mesh1 == mesh2);
}

BOOST_AUTO_TEST_CASE(MeshWKTPrint)
{
  PRECICE_TEST(1_rank);
  Mesh    mesh("WKTMesh", 3, testing::nextMeshID());
  Vertex &v0 = mesh.createVertex(Eigen::Vector3d(0., 0., 0.));
  Vertex &v1 = mesh.createVertex(Eigen::Vector3d(1., 0., 0.));
  Vertex &v2 = mesh.createVertex(Eigen::Vector3d(0., 0., 1.));
  Vertex &v3 = mesh.createVertex(Eigen::Vector3d(1., 0., 1.));
  Edge &  e0 = mesh.createEdge(v0, v1); // LINESTRING (0 0 0, 1 0 0)
  Edge &  e1 = mesh.createEdge(v1, v2); // LINESTRING (1 0 0, 0 0 1)
  Edge &  e2 = mesh.createEdge(v2, v0); // LINESTRING (0 0 1, 0 0 0)
  mesh.createEdge(v1, v3);              // LINESTRING (1 0 0, 1 0 1)
  mesh.createEdge(v3, v2);              // LINESTRING (1 0 1, 0 0 1)
  mesh.createTriangle(e0, e1, e2);
  std::stringstream sstream;
  sstream << mesh;
  std::string reference(
      "Mesh \"WKTMesh\", dimensionality = 3:\n"
      "GEOMETRYCOLLECTION(\n"
      "POINT (0 0 0), POINT (1 0 0), POINT (0 0 1), POINT (1 0 1),\n"
      "LINESTRING (0 0 0, 1 0 0), LINESTRING (1 0 0, 0 0 1), LINESTRING (0 0 1, 0 0 0), LINESTRING (1 0 0, 1 0 1), LINESTRING (1 0 1, 0 0 1),\n"
      "POLYGON ((0 0 0, 1 0 0, 0 0 1, 0 0 0))\n"
      ")");
  BOOST_TEST(reference == sstream.str());
}

BOOST_AUTO_TEST_CASE(CreateUniqueEdge)
{
  PRECICE_TEST(1_rank);
  int             dim = 3;
  Mesh            mesh1("Mesh1", dim, testing::nextMeshID());
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

BOOST_AUTO_TEST_CASE(ResizeDataGrow)
{
  PRECICE_TEST(1_rank);
  precice::mesh::Mesh mesh("MyMesh", 3, testing::nextMeshID());
  const auto &        values = mesh.createData("Data", 1, 0_dataID)->values();

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
  precice::mesh::Mesh mesh("MyMesh", 3, testing::nextMeshID());
  const auto &        values = mesh.createData("Data", 1, 0_dataID)->values();

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

BOOST_AUTO_TEST_SUITE(Utils)

BOOST_AUTO_TEST_CASE(AsChain)
{
  PRECICE_TEST(1_rank);
  Mesh mesh("Mesh1", 3, testing::nextMeshID());
  mesh.createData("Data", 1, 0_dataID);

  Eigen::Vector3d coords0;
  Eigen::Vector3d coords1;
  Eigen::Vector3d coords2;
  Eigen::Vector3d coords3;
  coords0 << 1.0, 0.0, 0.0;
  coords1 << 0.0, 0.0, 0.0;
  coords2 << 2.0, 2.0, 0.0;
  coords3 << 0.0, 1.0, 0.0;
  Vertex &v0 = mesh.createVertex(coords0);
  Vertex &v1 = mesh.createVertex(coords1);
  Vertex &v2 = mesh.createVertex(coords2);
  Vertex &v3 = mesh.createVertex(coords3);

  BOOST_TEST(mesh.vertices().size() == 4);

  Edge &e0 = mesh.createEdge(v0, v1);
  Edge &e1 = mesh.createEdge(v3, v2);
  Edge &e2 = mesh.createEdge(v3, v1);
  Edge &e3 = mesh.createEdge(v0, v2);

  BOOST_TEST(mesh.edges().size() == 4);

  auto chain = asChain(utils::make_array(&e0, &e1, &e2, &e3));

  BOOST_REQUIRE(chain.connected);

  BOOST_TEST(chain.edges.at(0) == &e0);
  BOOST_TEST(chain.edges.at(1) == &e2);
  BOOST_TEST(chain.edges.at(2) == &e1);
  BOOST_TEST(chain.edges.at(3) == &e3);

  BOOST_TEST(chain.vertices.at(0) == &v1);
  BOOST_TEST(chain.vertices.at(1) == &v3);
  BOOST_TEST(chain.vertices.at(2) == &v2);
  BOOST_TEST(chain.vertices.at(3) == &v0);
}

BOOST_AUTO_TEST_CASE(ShareVertex)
{
  PRECICE_TEST(1_rank);
  Mesh mesh("Mesh1", 3, testing::nextMeshID());
  mesh.createData("Data", 1, 0_dataID);

  Eigen::Vector3d coords0;
  Eigen::Vector3d coords1;
  Eigen::Vector3d coords2;
  Eigen::Vector3d coords3;
  coords0 << 1.0, 0.0, 0.0;
  coords1 << 0.0, 0.0, 0.0;
  coords2 << 2.0, 2.0, 0.0;
  coords3 << 0.0, 1.0, 0.0;
  Vertex *v0 = &mesh.createVertex(coords0);
  Vertex *v1 = &mesh.createVertex(coords1);
  Vertex *v2 = &mesh.createVertex(coords2);
  Vertex *v3 = &mesh.createVertex(coords3);
  BOOST_REQUIRE(mesh.vertices().size() == 4);

  Edge &e0 = mesh.createEdge(*v0, *v1);
  Edge &e1 = mesh.createEdge(*v1, *v2);
  Edge &e2 = mesh.createEdge(*v2, *v3);
  Edge &e3 = mesh.createEdge(*v3, *v0);
  BOOST_REQUIRE(mesh.edges().size() == 4);

  BOOST_TEST(sharedVertex(e0, e1) == v1);
  BOOST_TEST(sharedVertex(e1, e2) == v2);
  BOOST_TEST(sharedVertex(e2, e3) == v3);
  BOOST_TEST(sharedVertex(e3, e0) == v0);

  // not connnected
  BOOST_TEST((sharedVertex(e0, e2) == nullptr));
  BOOST_TEST((sharedVertex(e1, e3) == nullptr));

  // check symmetry
  BOOST_TEST(sharedVertex(e0, e1) == sharedVertex(e1, e0));
  BOOST_TEST(sharedVertex(e0, e2) == sharedVertex(e2, e0));
}

BOOST_AUTO_TEST_CASE(EdgeLength)
{
  PRECICE_TEST(1_rank);

  Eigen::Vector3d coords0;
  Eigen::Vector3d coords1;
  coords0 << 1.0, 0.0, 0.0;
  coords1 << 0.0, 1.0, 0.0;
  Vertex v0{coords0, 0};
  Vertex v1{coords1, 1};
  Edge   e(v0, v1, 0);
  BOOST_TEST(edgeLength(e) == std::sqrt(2));
}

BOOST_AUTO_TEST_CASE(VertexPtrsFor)
{
  PRECICE_TEST(1_rank);
  Mesh mesh("Mesh1", 3, testing::nextMeshID());
  mesh.createData("Data", 1, 0_dataID);

  Eigen::Vector3d coords0;
  Eigen::Vector3d coords1;
  Eigen::Vector3d coords2;
  Eigen::Vector3d coords3;
  coords0 << 1.0, 0.0, 0.0;
  coords1 << 0.0, 0.0, 0.0;
  coords2 << 2.0, 2.0, 0.0;
  coords3 << 0.0, 1.0, 0.0;
  Vertex &v0 = mesh.createVertex(coords0);
  Vertex &v1 = mesh.createVertex(coords1);
  Vertex &v2 = mesh.createVertex(coords2);
  mesh.createVertex(coords3);
  BOOST_TEST(mesh.vertices().size() == 4);

  std::array<int, 3>      ids{v0.getID(), v2.getID(), v1.getID()};
  std::array<Vertex *, 3> expected{&v0, &v2, &v1};

  auto result = vertexPtrsFor(mesh, ids);
  BOOST_REQUIRE(result.size() == 3);
  BOOST_TEST(result == expected);
}

BOOST_AUTO_TEST_CASE(CoordsForIDs)
{
  PRECICE_TEST(1_rank);
  Mesh mesh("Mesh1", 3, testing::nextMeshID());
  mesh.createData("Data", 1, 0_dataID);

  Eigen::Vector3d coords0;
  Eigen::Vector3d coords1;
  Eigen::Vector3d coords2;
  Eigen::Vector3d coords3;
  coords0 << 1.0, 0.0, 0.0;
  coords1 << 0.0, 0.0, 0.0;
  coords2 << 2.0, 2.0, 0.0;
  coords3 << 0.0, 1.0, 0.0;
  Vertex &v0 = mesh.createVertex(coords0);
  Vertex &v1 = mesh.createVertex(coords1);
  Vertex &v2 = mesh.createVertex(coords2);
  mesh.createVertex(coords3);
  BOOST_TEST(mesh.vertices().size() == 4);

  std::array<int, 3>             ids{v0.getID(), v2.getID(), v1.getID()};
  std::array<Eigen::VectorXd, 3> expected{coords0, coords2, coords1};

  auto result = coordsFor(mesh, ids);
  BOOST_REQUIRE(result.size() == 3);
  BOOST_TEST(result == expected);
}

BOOST_AUTO_TEST_CASE(CoordsForPtrs)
{
  PRECICE_TEST(1_rank);
  Mesh mesh("Mesh1", 3, testing::nextMeshID());
  mesh.createData("Data", 1, 0_dataID);

  Eigen::Vector3d coords0;
  Eigen::Vector3d coords1;
  Eigen::Vector3d coords2;
  Eigen::Vector3d coords3;
  coords0 << 1.0, 0.0, 0.0;
  coords1 << 0.0, 0.0, 0.0;
  coords2 << 2.0, 2.0, 0.0;
  coords3 << 0.0, 1.0, 0.0;
  Vertex &v0 = mesh.createVertex(coords0);
  Vertex &v1 = mesh.createVertex(coords1);
  Vertex &v2 = mesh.createVertex(coords2);
  mesh.createVertex(coords3);
  BOOST_TEST(mesh.vertices().size() == 4);

  std::array<Vertex *, 3>        ptrs{&v0, &v2, &v1};
  std::array<Eigen::VectorXd, 3> expected{coords0, coords2, coords1};

  auto result = coordsFor(ptrs);
  BOOST_REQUIRE(result.size() == 3);
  BOOST_TEST(result == expected);
}

BOOST_AUTO_TEST_CASE(Integrate2DScalarData)
{
  PRECICE_TEST(1_rank);
  PtrMesh mesh = std::make_shared<Mesh>("Mesh1", 2, testing::nextMeshID());
  mesh->createData("Data", 1, 0_dataID);

  auto &v1 = mesh->createVertex(Eigen::Vector2d(0.0, 0.0));
  auto &v2 = mesh->createVertex(Eigen::Vector2d(0.0, 0.5));
  auto &v3 = mesh->createVertex(Eigen::Vector2d(0.0, 1.5));
  auto &v4 = mesh->createVertex(Eigen::Vector2d(0.0, 3.5));
  mesh->allocateDataValues();

  mesh->createEdge(v1, v2); // Length = 0.5
  mesh->createEdge(v2, v3); // Length = 1.0
  mesh->createEdge(v3, v4); // Length = 2.0

  mesh->data(0)->values()(0) = 1.0;
  mesh->data(0)->values()(1) = 3.0;
  mesh->data(0)->values()(2) = 5.0;
  mesh->data(0)->values()(3) = 7.0;

  auto   result   = mesh::integrate(mesh, mesh->data(0));
  double expected = 17.0;
  BOOST_REQUIRE(result.size() == 1);
  BOOST_TEST(result(0) == expected);
}

BOOST_AUTO_TEST_CASE(Integrate2DVectorData)
{
  PRECICE_TEST(1_rank);
  PtrMesh mesh = std::make_shared<Mesh>("Mesh1", 2, testing::nextMeshID());
  mesh->createData("Data", 2, 0_dataID);

  auto &v1 = mesh->createVertex(Eigen::Vector2d(0.0, 0.0));
  auto &v2 = mesh->createVertex(Eigen::Vector2d(0.0, 0.5));
  auto &v3 = mesh->createVertex(Eigen::Vector2d(0.0, 1.5));
  auto &v4 = mesh->createVertex(Eigen::Vector2d(0.0, 3.5));
  mesh->allocateDataValues();

  mesh->createEdge(v1, v2); // Length = 0.5
  mesh->createEdge(v2, v3); // Length = 1.0
  mesh->createEdge(v3, v4); // Length = 2.0

  mesh->data(0)->values()(0) = 1.0;
  mesh->data(0)->values()(1) = 2.0;
  mesh->data(0)->values()(2) = 3.0;
  mesh->data(0)->values()(3) = 4.0;
  mesh->data(0)->values()(4) = 5.0;
  mesh->data(0)->values()(5) = 6.0;
  mesh->data(0)->values()(6) = 7.0;
  mesh->data(0)->values()(7) = 8.0;

  auto            result = mesh::integrate(mesh, mesh->data(0));
  Eigen::Vector2d expected(17.0, 20.5);
  BOOST_REQUIRE(result.size() == 2);
  BOOST_TEST(result(0) == expected(0));
  BOOST_TEST(result(1) == expected(1));
}

BOOST_AUTO_TEST_CASE(Integrate3DScalarData)
{
  PRECICE_TEST(1_rank);
  PtrMesh mesh = std::make_shared<Mesh>("Mesh1", 3, testing::nextMeshID());
  mesh->createData("Data", 1, 0_dataID);

  auto &v1 = mesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  auto &v2 = mesh->createVertex(Eigen::Vector3d(3.0, 0.0, 0.0));
  auto &v3 = mesh->createVertex(Eigen::Vector3d(3.0, 4.0, 0.0));
  auto &v4 = mesh->createVertex(Eigen::Vector3d(0.0, 8.0, 0.0));
  mesh->allocateDataValues();

  auto &e1 = mesh->createEdge(v1, v2);
  auto &e2 = mesh->createEdge(v2, v3);
  auto &e3 = mesh->createEdge(v3, v4);
  auto &e4 = mesh->createEdge(v1, v4);
  auto &e5 = mesh->createEdge(v1, v3);

  mesh->createTriangle(e1, e2, e5); // Area = 6.0
  mesh->createTriangle(e3, e4, e5); // Area = 12.0

  mesh->data(0)->values()(0) = 1.0;
  mesh->data(0)->values()(1) = 3.0;
  mesh->data(0)->values()(2) = 5.0;
  mesh->data(0)->values()(3) = 7.0;

  auto   result   = mesh::integrate(mesh, mesh->data(0));
  double expected = 70.0;
  BOOST_REQUIRE(result.size() == 1);
  BOOST_TEST(result(0) == expected);
}

BOOST_AUTO_TEST_CASE(Integrate3DVectorData)
{
  PRECICE_TEST(1_rank);
  PtrMesh mesh = std::make_shared<Mesh>("Mesh1", 3, testing::nextMeshID());
  mesh->createData("Data", 2, 0_dataID);

  auto &v1 = mesh->createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  auto &v2 = mesh->createVertex(Eigen::Vector3d(3.0, 0.0, 0.0));
  auto &v3 = mesh->createVertex(Eigen::Vector3d(3.0, 4.0, 0.0));
  auto &v4 = mesh->createVertex(Eigen::Vector3d(0.0, 8.0, 0.0));
  mesh->allocateDataValues();

  auto &e1 = mesh->createEdge(v1, v2);
  auto &e2 = mesh->createEdge(v2, v3);
  auto &e3 = mesh->createEdge(v3, v4);
  auto &e4 = mesh->createEdge(v1, v4);
  auto &e5 = mesh->createEdge(v1, v3);

  mesh->createTriangle(e1, e2, e5); // Area = 6.0
  mesh->createTriangle(e3, e4, e5); // Area = 12.0

  mesh->data(0)->values()(0) = 1.0;
  mesh->data(0)->values()(1) = 2.0;
  mesh->data(0)->values()(2) = 3.0;
  mesh->data(0)->values()(3) = 4.0;
  mesh->data(0)->values()(4) = 5.0;
  mesh->data(0)->values()(5) = 6.0;
  mesh->data(0)->values()(6) = 7.0;
  mesh->data(0)->values()(7) = 8.0;

  auto            result = mesh::integrate(mesh, mesh->data(0));
  Eigen::Vector2d expected(70.0, 88.0);
  BOOST_REQUIRE(result.size() == 2);
  BOOST_TEST(result(0) == expected(0));
  BOOST_TEST(result(1) == expected(1));
}

BOOST_AUTO_TEST_SUITE_END() // Utils

BOOST_AUTO_TEST_SUITE_END() // Mesh
BOOST_AUTO_TEST_SUITE_END() // Mesh
