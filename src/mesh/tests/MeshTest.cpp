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

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(BoundingBoxCOG_2D)
{
  PRECICE_TEST();
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

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(BoundingBoxCOG_3D)
{
  PRECICE_TEST();
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

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(ResetBoundingBox)
{
  PRECICE_TEST();
  // First test 2D mesh
  mesh::Mesh mesh2D("2D Testmesh", 2, testing::nextMeshID());
  {
    Eigen::Vector2d coords0(2, 0);
    Eigen::Vector2d coords1(-1, 4);
    Eigen::Vector2d coords2(0, 1);

    mesh2D.createVertex(coords0);
    mesh2D.createVertex(coords1);
    mesh2D.createVertex(coords2);
    mesh2D.computeBoundingBox();
  }

  const mesh::BoundingBox &bBox2D = mesh2D.getBoundingBox();
  mesh::BoundingBox        referenceBox2D({-1.0, 2.0,
                                    0.0, 4.0});

  BOOST_TEST(bBox2D.getDimension() == 2);
  BOOST_TEST(referenceBox2D == bBox2D);

  // The dimension remains
  mesh2D.resetBoundingBox();
  BOOST_TEST(bBox2D.getDimension() == 2);
  mesh::BoundingBox resetBox2D(bBox2D.getDimension());
  // Test that we reset the box
  BOOST_TEST(resetBox2D == bBox2D);

  // Now the 3D case
  mesh::Mesh mesh3D("3D Testmesh", 3, testing::nextMeshID());
  {
    Eigen::Vector3d coords0(2, 0, -3);
    Eigen::Vector3d coords1(-1, 4, 8);
    Eigen::Vector3d coords2(0, 1, -2);
    Eigen::Vector3d coords3(3.5, 2, -2);
    mesh3D.createVertex(coords0);
    mesh3D.createVertex(coords1);
    mesh3D.createVertex(coords2);
    mesh3D.createVertex(coords3);
  }

  mesh3D.computeBoundingBox();
  mesh::BoundingBox bBox3D = mesh3D.getBoundingBox();
  mesh::BoundingBox referenceBox3D({-1.0, 3.5,
                                    0.0, 4.0,
                                    -3.0, 8.0});

  BOOST_TEST(bBox3D.getDimension() == 3);
  BOOST_TEST(referenceBox3D == bBox3D);
  mesh3D.resetBoundingBox();
  mesh::BoundingBox resetBox3D(bBox3D.getDimension());
  // Test that we reset the box
  BOOST_TEST(resetBox3D == bBox3D);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(Demonstration)
{
  PRECICE_TEST();
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

    BOOST_TEST(!mesh.hasEdges());
    BOOST_TEST(!mesh.hasTriangles());

    // Create mesh edges
    Edge &e0 = mesh.createEdge(v0, v1);
    Edge &e1 = mesh.createEdge(v1, v2);
    Edge &e2 = mesh.createEdge(v2, v0);

    BOOST_TEST(mesh.hasEdges());
    BOOST_TEST(mesh.edges().size() == 3);
    BOOST_TEST(!mesh.hasTriangles());

    Triangle *t = nullptr;
    if (dim == 3) {
      // Create triangle
      t = &mesh.createTriangle(e0, e1, e2);

      BOOST_TEST(mesh.hasTriangles());
    } else {
      BOOST_TEST(!mesh.hasTriangles());
    }

    BOOST_TEST(mesh.hasEdges());

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

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(MeshEquality)
{
  PRECICE_TEST();
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

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(MeshWKTPrint)
{
  PRECICE_TEST();
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
      "LINESTRING (0 0 0, 1 0 0), LINESTRING (1 0 0, 0 0 1), LINESTRING (0 0 0, 0 0 1), LINESTRING (1 0 0, 1 0 1), LINESTRING (0 0 1, 1 0 1),\n"
      "POLYGON ((0 0 0, 1 0 0, 0 0 1, 0 0 0))\n"
      ")");
  BOOST_TEST(reference == sstream.str());
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(ResizeDataGrow)
{
  PRECICE_TEST();
  precice::mesh::Mesh mesh("MyMesh", 3, testing::nextMeshID());
  const auto &        values = mesh.createData("Data", 1, 0_dataID)->values();

  // Create mesh
  mesh.createVertex(Vector3d(0.0, 0.0, 0.0));
  mesh.createVertex(Vector3d(1.0, 0.0, 1.0));

  BOOST_TEST(mesh.nVertices() == 2);
  mesh.allocateDataValues();
  BOOST_TEST(values.size() == 2);

  mesh.createVertex(Vector3d(1.0, 1.0, 1.0));
  mesh.createVertex(Vector3d(2.0, 0.0, 2.0));
  mesh.createVertex(Vector3d(2.0, 0.0, 2.1));

  BOOST_TEST(mesh.nVertices() == 5);
  mesh.allocateDataValues();
  BOOST_TEST(values.size() == 5);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(ResizeDataShrink)
{
  PRECICE_TEST();
  precice::mesh::Mesh mesh("MyMesh", 3, testing::nextMeshID());
  const auto &        values = mesh.createData("Data", 1, 0_dataID)->values();

  // Create mesh
  mesh.createVertex(Vector3d(0.0, 0.0, 0.0));
  mesh.createVertex(Vector3d(1.0, 0.0, 1.0));
  mesh.createVertex(Vector3d(1.0, 1.0, 1.0));
  mesh.createVertex(Vector3d(2.0, 2.0, 2.0));

  BOOST_TEST(mesh.nVertices() == 4);
  mesh.allocateDataValues();
  BOOST_TEST(values.size() == 4);

  mesh.clear();
  mesh.createVertex(Vector3d(0.0, 0.0, 0.0));
  mesh.createVertex(Vector3d(1.0, 0.0, 1.0));

  BOOST_TEST(mesh.nVertices() == 2);
  mesh.allocateDataValues();
  BOOST_TEST(values.size() == 2);
}

BOOST_AUTO_TEST_SUITE(Utils)

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(AsChain)
{
  PRECICE_TEST();
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

  BOOST_TEST(mesh.nVertices() == 4);

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

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(ShareVertex)
{
  PRECICE_TEST();
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
  BOOST_REQUIRE(mesh.nVertices() == 4);

  Edge &e0 = mesh.createEdge(*v0, *v1);
  Edge &e1 = mesh.createEdge(*v1, *v2);
  Edge &e2 = mesh.createEdge(*v2, *v3);
  Edge &e3 = mesh.createEdge(*v3, *v0);
  BOOST_REQUIRE(mesh.edges().size() == 4);

  BOOST_TEST(sharedVertex(e0, e1) == v1);
  BOOST_TEST(sharedVertex(e1, e2) == v2);
  BOOST_TEST(sharedVertex(e2, e3) == v3);
  BOOST_TEST(sharedVertex(e3, e0) == v0);

  // not connected
  BOOST_TEST((sharedVertex(e0, e2) == nullptr));
  BOOST_TEST((sharedVertex(e1, e3) == nullptr));

  // check symmetry
  BOOST_TEST(sharedVertex(e0, e1) == sharedVertex(e1, e0));
  BOOST_TEST(sharedVertex(e0, e2) == sharedVertex(e2, e0));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(EdgeLength)
{
  PRECICE_TEST();

  Eigen::Vector3d coords0;
  Eigen::Vector3d coords1;
  coords0 << 1.0, 0.0, 0.0;
  coords1 << 0.0, 1.0, 0.0;
  Vertex v0{coords0, 0};
  Vertex v1{coords1, 1};
  Edge   e(v0, v1);
  BOOST_TEST(edgeLength(e) == std::sqrt(2));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(VertexPtrsFor)
{
  PRECICE_TEST();
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
  BOOST_TEST(mesh.nVertices() == 4);

  std::array<int, 3>      ids{v0.getID(), v2.getID(), v1.getID()};
  std::array<Vertex *, 3> expected{&v0, &v2, &v1};

  auto result = vertexPtrsFor(mesh, ids);
  BOOST_REQUIRE(result.size() == 3);
  BOOST_TEST(result == expected);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(CoordsForIDs)
{
  PRECICE_TEST();
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
  BOOST_TEST(mesh.nVertices() == 4);

  std::array<int, 3>             ids{v0.getID(), v2.getID(), v1.getID()};
  std::array<Eigen::VectorXd, 3> expected{coords0, coords2, coords1};

  auto result = coordsFor(mesh, ids);
  BOOST_REQUIRE(result.size() == 3);
  BOOST_TEST(result == expected);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(CoordsForPtrs)
{
  PRECICE_TEST();
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
  BOOST_TEST(mesh.nVertices() == 4);

  std::array<Vertex *, 3>        ptrs{&v0, &v2, &v1};
  std::array<Eigen::VectorXd, 3> expected{coords0, coords2, coords1};

  auto result = coordsFor(ptrs);
  BOOST_REQUIRE(result.size() == 3);
  BOOST_TEST(result == expected);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(Integrate2DScalarData)
{
  PRECICE_TEST();
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

  auto   result   = mesh::integrateSurface(mesh, mesh->data(0)->values());
  double expected = 17.0;
  BOOST_REQUIRE(result.size() == 1);
  BOOST_TEST(result(0) == expected);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(Integrate2DVectorData)
{
  PRECICE_TEST();
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

  auto            result = mesh::integrateSurface(mesh, mesh->data(0)->values());
  Eigen::Vector2d expected(17.0, 20.5);
  BOOST_REQUIRE(result.size() == 2);
  BOOST_TEST(result(0) == expected(0));
  BOOST_TEST(result(1) == expected(1));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(Integrate3DScalarData)
{
  PRECICE_TEST();
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

  auto   result   = mesh::integrateSurface(mesh, mesh->data(0)->values());
  double expected = 70.0;
  BOOST_REQUIRE(result.size() == 1);
  BOOST_TEST(result(0) == expected);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(Integrate3DVectorData)
{
  PRECICE_TEST();
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

  auto            result = mesh::integrateSurface(mesh, mesh->data(0)->values());
  Eigen::Vector2d expected(70.0, 88.0);
  BOOST_REQUIRE(result.size() == 2);
  BOOST_TEST(result(0) == expected(0));
  BOOST_TEST(result(1) == expected(1));
}

BOOST_AUTO_TEST_SUITE_END() // Utils

BOOST_AUTO_TEST_SUITE(VolumeIntegrals)

struct UnitSquareFixture {
  Eigen::Vector2d x0{0.0, 0.0};
  Eigen::Vector2d x1{1.0, 0.0};
  Eigen::Vector2d x2{1.0, 1.0};
  Eigen::Vector2d x3{0.0, 1.0};
};

PRECICE_TEST_SETUP(1_rank)
BOOST_FIXTURE_TEST_CASE(Integrate2DScalarDataVolume, UnitSquareFixture)
{
  PRECICE_TEST();
  PtrMesh mesh = std::make_shared<Mesh>("Mesh1", 2, testing::nextMeshID());
  mesh->createData("Data", 1, 0_dataID);

  auto &v1 = mesh->createVertex(x0);
  auto &v2 = mesh->createVertex(x1);
  auto &v3 = mesh->createVertex(x2);
  auto &v4 = mesh->createVertex(x3);
  mesh->allocateDataValues();

  auto &e1 = mesh->createEdge(v1, v2);
  auto &e2 = mesh->createEdge(v2, v3);
  auto &e3 = mesh->createEdge(v3, v4);
  auto &e4 = mesh->createEdge(v1, v4);
  auto &e5 = mesh->createEdge(v1, v3);

  // 2 triangles as halves of a unit square
  mesh->createTriangle(e1, e2, e5);
  mesh->createTriangle(e3, e4, e5);

  BOOST_REQUIRE(mesh->triangles()[0].getArea() == 0.5);
  BOOST_REQUIRE(mesh->triangles()[1].getArea() == 0.5);
  BOOST_REQUIRE(mesh->triangles().size() == 2);

  // Integrand is 1 + 4x + 2y integrated over unit square: 4
  mesh->data(0)->values()(0) = 1.0;
  mesh->data(0)->values()(1) = 3.0;
  mesh->data(0)->values()(2) = 7.0;
  mesh->data(0)->values()(3) = 5.0;

  auto   result   = mesh::integrateVolume(mesh, mesh->data(0)->values());
  double expected = 4.0;
  BOOST_REQUIRE(result.size() == 1);
  BOOST_TEST(result(0) == expected);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_FIXTURE_TEST_CASE(Integrate2DVectorDataVolume, UnitSquareFixture)
{
  PRECICE_TEST();
  PtrMesh mesh = std::make_shared<Mesh>("Mesh1", 2, testing::nextMeshID());
  mesh->createData("Data", 2, 0_dataID);

  auto &v1 = mesh->createVertex(x0);
  auto &v2 = mesh->createVertex(x1);
  auto &v3 = mesh->createVertex(x2);
  auto &v4 = mesh->createVertex(x3);
  mesh->allocateDataValues();

  auto &e1 = mesh->createEdge(v1, v2);
  auto &e2 = mesh->createEdge(v2, v3);
  auto &e3 = mesh->createEdge(v3, v4);
  auto &e4 = mesh->createEdge(v1, v4);
  auto &e5 = mesh->createEdge(v1, v3);

  // 2 triangles as halves of a unit square
  mesh->createTriangle(e1, e2, e5);
  mesh->createTriangle(e3, e4, e5);

  BOOST_REQUIRE(mesh->triangles()[0].getArea() == 0.5);
  BOOST_REQUIRE(mesh->triangles()[1].getArea() == 0.5);
  BOOST_REQUIRE(mesh->triangles().size() == 2);

  // Integrand of 1st component is 1 + 4x + 2y integrated over unit square: 4
  mesh->data(0)->values()(0) = 1.0;
  mesh->data(0)->values()(2) = 3.0;
  mesh->data(0)->values()(4) = 7.0;
  mesh->data(0)->values()(6) = 5.0;
  // Integrand of 1st component is 1 + 4x + 2y integrated over unit square: 4
  mesh->data(0)->values()(1) = 2.0;
  mesh->data(0)->values()(3) = 4.0;
  mesh->data(0)->values()(5) = 8.0;
  mesh->data(0)->values()(7) = 6.0;

  auto            result = mesh::integrateVolume(mesh, mesh->data(0)->values());
  Eigen::Vector2d expected(4.0, 5.0);
  BOOST_REQUIRE(result.size() == 2);
  BOOST_TEST(result(0) == expected(0));
  BOOST_TEST(result(1) == expected(1));
}

struct OneTetraFixture {
  Eigen::Vector3d x1{0.0, 0.0, 0.0};
  Eigen::Vector3d x2{1.0, 0.0, 0.0};
  Eigen::Vector3d x3{0.0, 1.0, 0.0};
  Eigen::Vector3d x4{0.0, 0.0, 1.0};
};

PRECICE_TEST_SETUP(1_rank)
BOOST_FIXTURE_TEST_CASE(Integrate3DScalarDataVolume, OneTetraFixture)
{
  PRECICE_TEST();
  PtrMesh mesh = std::make_shared<Mesh>("Mesh1", 3, testing::nextMeshID());
  mesh->createData("Data", 1, 0_dataID);

  auto &v1 = mesh->createVertex(x1);
  auto &v2 = mesh->createVertex(x2);
  auto &v3 = mesh->createVertex(x3);
  auto &v4 = mesh->createVertex(x4);

  mesh->allocateDataValues();

  mesh->createTetrahedron(v1, v2, v3, v4);

  BOOST_REQUIRE(mesh->tetrahedra()[0].getVolume() == 1. / 6);
  BOOST_REQUIRE(mesh->tetrahedra().size() == 1);

  // Integrand is 1 + 2x + 4y + 6z integrated over one tetra
  mesh->data(0)->values()(0) = 1.0;
  mesh->data(0)->values()(1) = 3.0;
  mesh->data(0)->values()(2) = 5.0;
  mesh->data(0)->values()(3) = 7.0;

  auto   result   = mesh::integrateVolume(mesh, mesh->data(0)->values());
  double expected = 4.0 / 6;
  BOOST_REQUIRE(result.size() == 1);
  BOOST_TEST(result(0) == expected);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_FIXTURE_TEST_CASE(Integrate3DVectorDataVolume, OneTetraFixture)
{
  PRECICE_TEST();
  PtrMesh mesh = std::make_shared<Mesh>("Mesh1", 3, testing::nextMeshID());
  mesh->createData("Data", 3, 0_dataID);

  auto &v1 = mesh->createVertex(x1);
  auto &v2 = mesh->createVertex(x2);
  auto &v3 = mesh->createVertex(x3);
  auto &v4 = mesh->createVertex(x4);

  mesh->allocateDataValues();

  mesh->createTetrahedron(v1, v2, v3, v4);

  BOOST_REQUIRE(mesh->tetrahedra()[0].getVolume() == 1. / 6);
  BOOST_REQUIRE(mesh->tetrahedra().size() == 1);

  // Integrand is (1 + 2x + 4y + 6z, 1, 1-x-y-z) integrated over one tetra
  mesh->data(0)->values()(0 * 3) = 1.0;
  mesh->data(0)->values()(1 * 3) = 3.0;
  mesh->data(0)->values()(2 * 3) = 5.0;
  mesh->data(0)->values()(3 * 3) = 7.0;

  mesh->data(0)->values()(0 * 3 + 1) = 1.0;
  mesh->data(0)->values()(1 * 3 + 1) = 1.0;
  mesh->data(0)->values()(2 * 3 + 1) = 1.0;
  mesh->data(0)->values()(3 * 3 + 1) = 1.0;

  mesh->data(0)->values()(0 * 3 + 2) = 1.0;
  mesh->data(0)->values()(1 * 3 + 2) = 0.0;
  mesh->data(0)->values()(2 * 3 + 2) = 0.0;
  mesh->data(0)->values()(3 * 3 + 2) = 0.0;

  auto            result = mesh::integrateVolume(mesh, mesh->data(0)->values());
  Eigen::Vector3d expected(4.0 / 6, 1.0 / 6, 1.0 / 24);
  BOOST_REQUIRE(result.size() == 3);
  BOOST_TEST(result(0) == expected(0));
  BOOST_TEST(result(1) == expected(1));
  BOOST_TEST(result(2) == expected(2));
}

BOOST_AUTO_TEST_SUITE_END() // VolumeIntegrals

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(AddMesh)
{
  PRECICE_TEST();
  PtrMesh globalMesh = std::make_shared<Mesh>("Mesh1", 3, testing::nextMeshID());
  PtrMesh subMesh    = std::make_shared<Mesh>("Mesh2", 3, testing::nextMeshID());

  // Fill globalMesh, then subMesh, then add and check the result is the sul
  auto &v01 = globalMesh->createVertex(Eigen::Vector3d{0.0, 0.0, 0.0});
  auto &v02 = globalMesh->createVertex(Eigen::Vector3d{1.0, 0.0, 0.0});
  auto &v03 = globalMesh->createVertex(Eigen::Vector3d{0.0, 1.0, 0.0});
  auto &v04 = globalMesh->createVertex(Eigen::Vector3d{0.0, 0.0, 1.0});

  // Add some elements of all kind to global
  globalMesh->createEdge(v01, v02);
  globalMesh->createTriangle(v01, v02, v03);
  globalMesh->createTetrahedron(v01, v02, v03, v04);

  auto &v11 = subMesh->createVertex(Eigen::Vector3d{0.0, 0.0, 0.0});
  auto &v12 = subMesh->createVertex(Eigen::Vector3d{2.0, 0.0, 0.0});
  auto &v13 = subMesh->createVertex(Eigen::Vector3d{0.0, 4.0, 0.0});
  auto &v14 = subMesh->createVertex(Eigen::Vector3d{0.0, 0.0, 3.0});
  auto &v15 = subMesh->createVertex(Eigen::Vector3d{5.0, 0.0, 1.0});

  // Add some elements of all kind to sub mesh
  subMesh->createEdge(v11, v12);
  subMesh->createEdge(v14, v15);
  subMesh->createTriangle(v11, v13, v15);
  subMesh->createTriangle(v11, v13, v14);
  subMesh->createTetrahedron(v11, v12, v13, v14);
  subMesh->createTetrahedron(v15, v12, v13, v14);

  globalMesh->addMesh(*subMesh);
  BOOST_TEST(globalMesh->nVertices() == 9);
  BOOST_TEST(globalMesh->edges().size() == 3);
  BOOST_TEST(globalMesh->triangles().size() == 3);
  BOOST_TEST(globalMesh->tetrahedra().size() == 3);
}

BOOST_AUTO_TEST_SUITE(PreProcess);

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(DuplicateEdges)
{
  PRECICE_TEST();
  Mesh mesh{"Mesh1", 3, 0};

  auto &v1 = mesh.createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  auto &v2 = mesh.createVertex(Eigen::Vector3d(3.0, 0.0, 0.0));
  auto &v3 = mesh.createVertex(Eigen::Vector3d(3.0, 4.0, 0.0));
  auto &v4 = mesh.createVertex(Eigen::Vector3d(0.0, 8.0, 0.0));

  mesh.createEdge(v1, v2);
  mesh.createEdge(v2, v3);
  mesh.createEdge(v3, v4);

  // add single duplicate
  mesh.createEdge(v3, v4);

  // add 1000 duplicates
  for (int i = 0; i < 1000; ++i) {
    mesh.createEdge(v2, v3);
  };
  BOOST_TEST(mesh.edges().size() == 1004);

  mesh.preprocess();

  BOOST_TEST(mesh.edges().size() == 3);
  BOOST_TEST(mesh.triangles().empty());
  BOOST_TEST(mesh.tetrahedra().empty());

  std::vector<Edge> expectedEdges{{v1, v2}, {v2, v3}, {v3, v4}};
  for (auto &e : expectedEdges) {
    auto cnt = std::count(mesh.edges().begin(), mesh.edges().end(), e);
    BOOST_TEST(cnt == 1);
  }
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(SingleTriangle)
{
  PRECICE_TEST();
  Mesh mesh{"Mesh1", 3, 0};

  auto &v1 = mesh.createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  auto &v2 = mesh.createVertex(Eigen::Vector3d(0.0, 1.0, 0.0));
  auto &v3 = mesh.createVertex(Eigen::Vector3d(1.0, 1.0, 0.0));

  mesh.createTriangle(v1, v2, v3);

  mesh.preprocess();

  BOOST_TEST(mesh.edges().size() == 3);
  BOOST_TEST(mesh.triangles().size() == 1);
  BOOST_TEST(mesh.tetrahedra().size() == 0);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(SingleTriangulatedQuad)
{
  PRECICE_TEST();
  Mesh mesh{"Mesh1", 3, 0};

  auto &v1 = mesh.createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  auto &v2 = mesh.createVertex(Eigen::Vector3d(1.0, 0.0, 0.0));
  auto &v3 = mesh.createVertex(Eigen::Vector3d(1.0, 1.0, 0.0));
  auto &v4 = mesh.createVertex(Eigen::Vector3d(0.0, 1.0, 0.0));

  mesh.createTriangle(v1, v2, v3);
  mesh.createTriangle(v2, v3, v4);

  mesh.preprocess();

  BOOST_TEST(mesh.edges().size() == 5);
  BOOST_TEST(mesh.triangles().size() == 2);
  BOOST_TEST(mesh.tetrahedra().size() == 0);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(SingleTetrahedron)
{
  PRECICE_TEST();
  Mesh mesh{"Mesh1", 3, 0};

  auto &v1 = mesh.createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  auto &v2 = mesh.createVertex(Eigen::Vector3d(0.0, 1.0, 0.0));
  auto &v3 = mesh.createVertex(Eigen::Vector3d(1.0, 1.0, 0.0));
  auto &v4 = mesh.createVertex(Eigen::Vector3d(0.3, 0.3, 1.0));

  mesh.createTetrahedron(v1, v2, v3, v4);

  mesh.preprocess();

  BOOST_TEST(mesh.edges().size() == 6);
  BOOST_TEST(mesh.triangles().size() == 4);
  BOOST_TEST(mesh.tetrahedra().size() == 1);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(TouchingTetrahedra)
{
  PRECICE_TEST();
  Mesh mesh{"Mesh1", 3, 0};

  auto &v1 = mesh.createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  auto &v2 = mesh.createVertex(Eigen::Vector3d(0.0, 1.0, 0.0));
  auto &v3 = mesh.createVertex(Eigen::Vector3d(1.0, 1.0, 0.0));

  auto &v4 = mesh.createVertex(Eigen::Vector3d(0.3, 0.3, 1.0));
  auto &v5 = mesh.createVertex(Eigen::Vector3d(0.3, 0.3, -1.0));

  mesh.createTetrahedron(v1, v2, v3, v4);
  mesh.createTetrahedron(v1, v2, v3, v5);

  mesh.preprocess();

  BOOST_TEST(mesh.edges().size() == 9);
  BOOST_TEST(mesh.triangles().size() == 7);
  BOOST_TEST(mesh.tetrahedra().size() == 2);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(SingleTriangle2ImplicitEdges)
{
  PRECICE_TEST();
  Mesh mesh{"Mesh1", 3, 0};

  auto &v1 = mesh.createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  auto &v2 = mesh.createVertex(Eigen::Vector3d(0.0, 1.0, 0.0));
  auto &v3 = mesh.createVertex(Eigen::Vector3d(1.0, 1.0, 0.0));

  mesh.createTriangle(v1, v2, v3);
  mesh.createEdge(v2, v3);

  mesh.preprocess();

  BOOST_TEST(mesh.edges().size() == 3);
  BOOST_TEST(mesh.triangles().size() == 1);
  BOOST_TEST(mesh.tetrahedra().size() == 0);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(SingleTriangle1ImplicitEdges)
{
  PRECICE_TEST();
  Mesh mesh{"Mesh1", 3, 0};

  auto &v1 = mesh.createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  auto &v2 = mesh.createVertex(Eigen::Vector3d(0.0, 1.0, 0.0));
  auto &v3 = mesh.createVertex(Eigen::Vector3d(1.0, 1.0, 0.0));

  mesh.createTriangle(v1, v2, v3);
  mesh.createEdge(v1, v2);
  mesh.createEdge(v2, v3);

  mesh.preprocess();

  BOOST_TEST(mesh.edges().size() == 3);
  BOOST_TEST(mesh.triangles().size() == 1);
  BOOST_TEST(mesh.tetrahedra().size() == 0);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(SingleTetrahedron3ImplicitTriangles)
{
  PRECICE_TEST();
  Mesh mesh{"Mesh1", 3, 0};

  auto &v1 = mesh.createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  auto &v2 = mesh.createVertex(Eigen::Vector3d(0.0, 1.0, 0.0));
  auto &v3 = mesh.createVertex(Eigen::Vector3d(1.0, 1.0, 0.0));
  auto &v4 = mesh.createVertex(Eigen::Vector3d(0.3, 0.3, 1.0));

  mesh.createTetrahedron(v1, v2, v3, v4);
  mesh.createTriangle(v1, v2, v3);

  mesh.preprocess();

  BOOST_TEST(mesh.edges().size() == 6);
  BOOST_TEST(mesh.triangles().size() == 4);
  BOOST_TEST(mesh.tetrahedra().size() == 1);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(SingleTetrahedronNoImplicitTriangles)
{
  PRECICE_TEST();
  Mesh mesh{"Mesh1", 3, 0};

  auto &v1 = mesh.createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  auto &v2 = mesh.createVertex(Eigen::Vector3d(0.0, 1.0, 0.0));
  auto &v3 = mesh.createVertex(Eigen::Vector3d(1.0, 1.0, 0.0));
  auto &v4 = mesh.createVertex(Eigen::Vector3d(0.3, 0.3, 1.0));

  mesh.createTetrahedron(v1, v2, v3, v4);
  mesh.createTriangle(v1, v2, v3);
  mesh.createTriangle(v1, v2, v4);
  mesh.createTriangle(v3, v4, v1);
  mesh.createTriangle(v3, v4, v2);

  mesh.preprocess();

  BOOST_TEST(mesh.edges().size() == 6);
  BOOST_TEST(mesh.triangles().size() == 4);
  BOOST_TEST(mesh.tetrahedra().size() == 1);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(Mixed)
{
  PRECICE_TEST();
  Mesh mesh{"Mesh1", 3, 0};

  auto &v1 = mesh.createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  auto &v2 = mesh.createVertex(Eigen::Vector3d(0.0, 1.0, 0.0));
  auto &v3 = mesh.createVertex(Eigen::Vector3d(1.0, 1.0, 0.0));
  auto &v4 = mesh.createVertex(Eigen::Vector3d(0.3, 0.3, 1.0));

  mesh.createTetrahedron(v1, v2, v3, v4);
  mesh.createTriangle(v1, v2, v3);
  mesh.createEdge(v3, v4);

  mesh.preprocess();

  BOOST_TEST(mesh.edges().size() == 6);
  BOOST_TEST(mesh.triangles().size() == 4);
  BOOST_TEST(mesh.tetrahedra().size() == 1);
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(Complex)
{
  PRECICE_TEST();
  Mesh mesh{"Mesh1", 3, 0};

  auto &v1 = mesh.createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  auto &v2 = mesh.createVertex(Eigen::Vector3d(3.0, 0.0, 0.0));
  auto &v3 = mesh.createVertex(Eigen::Vector3d(3.0, 4.0, 0.0));
  auto &v4 = mesh.createVertex(Eigen::Vector3d(0.0, 8.0, 0.0));

  mesh.createEdge(v1, v2);
  mesh.createEdge(v2, v3);
  mesh.createEdge(v3, v4);

  // add single duplicate
  mesh.createEdge(v3, v4);

  // add 1000 duplicates
  for (int i = 0; i < 1000; ++i) {
    mesh.createEdge(v2, v3);
  };
  BOOST_TEST(mesh.edges().size() == 1004);

  mesh.createTriangle(v1, v2, v3);
  mesh.createTriangle(v2, v3, v4);
  for (int i = 0; i < 1000; ++i) {
    mesh.createTriangle(v2, v3, v4);
  };
  BOOST_TEST(mesh.triangles().size() == 1002);

  // creates edges 1-3 and 2-4
  mesh.preprocess();

  BOOST_TEST(mesh.edges().size() == 5);
  BOOST_TEST(mesh.triangles().size() == 2);

  std::vector<Edge> expectedEdges{{v1, v2}, {v1, v3}, {v2, v3}, {v3, v4}, {v2, v4}};
  for (auto &e : expectedEdges) {
    auto cnt = std::count(mesh.edges().begin(), mesh.edges().end(), e);
    BOOST_TEST(cnt == 1);
  }

  std::vector<Triangle> expectedTriangles{{v1, v2, v3}, {v2, v3, v4}};
  for (auto &t : expectedTriangles) {
    auto cnt = std::count(mesh.triangles().begin(), mesh.triangles().end(), t);
    BOOST_TEST(cnt == 1);
  }
}

BOOST_AUTO_TEST_SUITE_END();
BOOST_AUTO_TEST_SUITE_END() // Mesh
BOOST_AUTO_TEST_SUITE_END() // Mesh
