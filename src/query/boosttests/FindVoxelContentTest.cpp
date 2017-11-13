#include "io/ExportVTK.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Vertex.hpp"
#include "query/ExportVTKVoxelQueries.hpp"
#include "query/FindVoxelContent.hpp"
#include "testing/Testing.hpp"
#include "utils/Globals.hpp"

using namespace precice;
using namespace precice::query;

BOOST_AUTO_TEST_SUITE(QueryTests)
BOOST_AUTO_TEST_SUITE(FindVoxelContentTests)

void performTestVertices(
    int                    testDim,
    bool                   positive,
    const Eigen::VectorXd &offset)
{
  int dim = offset.size();
  assertion(not math::oneGreater(offset, Eigen::VectorXd::Constant(dim, 1.0)));
  assertion(math::allGreater(offset, Eigen::VectorXd::Constant(dim, -1.0)));
  bool            flipNormals = false;
  mesh::Mesh      mesh("TestMesh", dim, flipNormals);
  Eigen::VectorXd coords(offset);
  mesh::Vertex &  vertex = mesh.createVertex(coords);

  Eigen::VectorXd                     center        = Eigen::VectorXd::Zero(dim);
  Eigen::VectorXd                     halflengths   = Eigen::VectorXd::Constant(dim, 1.0);
  FindVoxelContent::BoundaryInclusion includeBounds = FindVoxelContent::INCLUDE_BOUNDARY;
  FindVoxelContent::BoundaryInclusion excludeBounds = FindVoxelContent::EXCLUDE_BOUNDARY;
  query::FindVoxelContent             findIncluded(center, halflengths, includeBounds);
  query::FindVoxelContent             findExcluded(center, halflengths, excludeBounds);

  assertion(testDim >= 0);
  assertion(testDim < dim);

  double sign = positive ? 1.0 : -1.0;
  int    size = 0;

  // Outside
  coords[testDim] = sign * 2.0;
  vertex.setCoords(coords);
  findIncluded(mesh);
  findExcluded(mesh);
  size = findIncluded.content().vertices().size();
  BOOST_TEST(size == 0);
  size = findExcluded.content().vertices().size();
  BOOST_TEST(size == 0);

  // Outside eps
  coords[testDim] = sign * (1.0 + 10.0 * math::NUMERICAL_ZERO_DIFFERENCE);
  vertex.setCoords(coords);
  findIncluded(mesh);
  findExcluded(mesh);
  size = findIncluded.content().vertices().size();
  BOOST_TEST(size == 0);
  size = findExcluded.content().vertices().size();
  BOOST_TEST(size == 0);

  // Outside eps
  coords[testDim] = sign * (1.0 + 10.0 * math::NUMERICAL_ZERO_DIFFERENCE);
  vertex.setCoords(coords);
  findIncluded(mesh);
  findExcluded(mesh);
  size = findIncluded.content().vertices().size();
  BOOST_TEST(size == 0);
  size = findExcluded.content().vertices().size();
  BOOST_TEST(size == 0);

  // Touching + eps
  coords[testDim] = sign * (1.0 + math::NUMERICAL_ZERO_DIFFERENCE);
  vertex.setCoords(coords);
  findIncluded(mesh);
  findExcluded(mesh);
  size = findIncluded.content().vertices().size();
  BOOST_TEST(size == 1);
  size = findExcluded.content().vertices().size();
  BOOST_TEST(size == 0);

  // Touching
  coords[testDim] = sign * 1.0;
  vertex.setCoords(coords);
  findIncluded(mesh);
  findExcluded(mesh);
  size = findIncluded.content().vertices().size();
  BOOST_TEST(size == 2);
  size = findExcluded.content().vertices().size();
  BOOST_TEST(size == 0);

  // Touching - eps
  coords[testDim] = sign * (1.0 - math::NUMERICAL_ZERO_DIFFERENCE);
  vertex.setCoords(coords);
  findIncluded(mesh);
  findExcluded(mesh);
  size = findIncluded.content().vertices().size();
  BOOST_TEST(size == 3);
  size = findExcluded.content().vertices().size();
  BOOST_TEST(size == 0);

  // Inside eps
  coords[testDim] = sign * (1.0 - 10.0 * math::NUMERICAL_ZERO_DIFFERENCE);
  vertex.setCoords(coords);
  findIncluded(mesh);
  findExcluded(mesh);
  size = findIncluded.content().vertices().size();
  BOOST_TEST(size == 4);
  size = findExcluded.content().vertices().size();
  BOOST_TEST(size == 1);

  // Inside
  coords[testDim] = sign * 0.9;
  vertex.setCoords(coords);
  findIncluded(mesh);
  findExcluded(mesh);
  size = findIncluded.content().vertices().size();
  BOOST_TEST(size == 5);
  size = findExcluded.content().vertices().size();
  BOOST_TEST(size == 2);
}

void performTestEdges(
    int                    testDim,
    bool                   positive,
    const Eigen::VectorXd &offset)
{
  int dim = offset.size();
  assertion(not math::oneGreater(offset, Eigen::VectorXd::Constant(dim, 1.0)));
  assertion(math::allGreater(offset, Eigen::VectorXd::Constant(dim, -1.0)));
  bool            flipNormals = false;
  mesh::Mesh      mesh("TestMesh", dim, flipNormals);
  Eigen::VectorXd coords0(offset);
  Eigen::VectorXd coords1(offset);
  mesh::Vertex &  v0 = mesh.createVertex(coords0);
  mesh::Vertex &  v1 = mesh.createVertex(coords1);
  mesh.createEdge(v0, v1);

  Eigen::VectorXd                     center        = Eigen::VectorXd::Zero(dim);
  Eigen::VectorXd                     halflengths   = Eigen::VectorXd::Constant(dim, 1);
  FindVoxelContent::BoundaryInclusion includeBounds = FindVoxelContent::INCLUDE_BOUNDARY;
  FindVoxelContent::BoundaryInclusion excludeBounds = FindVoxelContent::EXCLUDE_BOUNDARY;
  query::FindVoxelContent             findIncluded(center, halflengths, includeBounds);
  query::FindVoxelContent             findExcluded(center, halflengths, excludeBounds);

  assertion(testDim >= 0);
  assertion(testDim < dim);

  double sign = positive ? 1.0 : -1.0;

  // Outside
  coords0[testDim] = sign * 2.0;
  coords1[testDim] = sign * 3.0;
  v0.setCoords(coords0);
  v1.setCoords(coords1);
  mesh.computeState();
  findIncluded(mesh);
  findExcluded(mesh);
  int sizeEdges = findIncluded.content().edges().size();
  BOOST_TEST(sizeEdges == 0);
  sizeEdges = findExcluded.content().edges().size();
  BOOST_TEST(sizeEdges == 0);

  // Outside eps
  coords0[testDim] = sign * (1.0 + 10.0 * math::NUMERICAL_ZERO_DIFFERENCE);
  coords1[testDim] = sign * 2.0;
  v0.setCoords(coords0);
  v1.setCoords(coords1);
  mesh.computeState();
  findIncluded(mesh);
  findExcluded(mesh);
  sizeEdges = findIncluded.content().edges().size();
  BOOST_TEST(sizeEdges == 0);
  sizeEdges = findExcluded.content().edges().size();
  BOOST_TEST(sizeEdges == 0);

  // Outside touching
  coords0[testDim] = sign * 1.0;
  coords1[testDim] = sign * 2.0;
  v0.setCoords(coords0);
  v1.setCoords(coords1);
  mesh.computeState();
  findIncluded(mesh);
  findExcluded(mesh);
  sizeEdges = findIncluded.content().edges().size();
  BOOST_TEST(sizeEdges == 1);
  sizeEdges = findExcluded.content().edges().size();
  BOOST_TEST(sizeEdges == 0);

  // Outside touching eps
  coords0[testDim] = sign * (1.0 + math::NUMERICAL_ZERO_DIFFERENCE);
  coords1[testDim] = sign * 2.0;
  v0.setCoords(coords0);
  v1.setCoords(coords1);
  mesh.computeState();
  findIncluded(mesh);
  findExcluded(mesh);
  sizeEdges = findIncluded.content().edges().size();
  BOOST_TEST(sizeEdges == 2);
  sizeEdges = findExcluded.content().edges().size();
  BOOST_TEST(sizeEdges == 0);

  // Intersecting
  coords0[testDim] = sign * 0.5;
  coords1[testDim] = sign * 1.5;
  v0.setCoords(coords0);
  v1.setCoords(coords1);
  mesh.computeState();
  findIncluded(mesh);
  findExcluded(mesh);
  sizeEdges = findIncluded.content().edges().size();
  BOOST_TEST(sizeEdges == 3);
  sizeEdges = findExcluded.content().edges().size();
  BOOST_TEST(sizeEdges == 1);

  // Intersecting eps
  coords0[testDim] = sign * (1.0 - 10.0 * math::NUMERICAL_ZERO_DIFFERENCE);
  coords1[testDim] = sign * 2.0;
  v0.setCoords(coords0);
  v1.setCoords(coords1);
  mesh.computeState();
  findIncluded(mesh);
  findExcluded(mesh);
  sizeEdges = findIncluded.content().edges().size();
  BOOST_TEST(sizeEdges == 4);
  sizeEdges = findExcluded.content().edges().size();
  BOOST_TEST(sizeEdges == 2);

  // Inside
  coords0[testDim] = sign * 0.3;
  coords1[testDim] = sign * 0.7;
  v0.setCoords(coords0);
  v1.setCoords(coords1);
  mesh.computeState();
  findIncluded(mesh);
  findExcluded(mesh);
  sizeEdges = findIncluded.content().edges().size();
  BOOST_TEST(sizeEdges == 5);
  sizeEdges = findExcluded.content().edges().size();
  BOOST_TEST(sizeEdges == 3);
}

void performTestTriangles(
    int  testDim,
    int  secondDimension,
    bool positive)
{
  int dim = 3;
  assertion(testDim != secondDimension);
  bool            flipNormals = false;
  mesh::Mesh      mesh("TestMesh", dim, flipNormals);
  Eigen::Vector3d coords0 = Eigen::Vector3d::Zero();
  Eigen::Vector3d coords1 = Eigen::Vector3d::Zero();
  Eigen::Vector3d coords2 = Eigen::Vector3d::Zero();
  mesh::Vertex &  v0      = mesh.createVertex(coords0);
  mesh::Vertex &  v1      = mesh.createVertex(coords1);
  mesh::Vertex &  v2      = mesh.createVertex(coords2);
  mesh::Edge &    e0      = mesh.createEdge(v0, v1);
  mesh::Edge &    e1      = mesh.createEdge(v1, v2);
  mesh::Edge &    e2      = mesh.createEdge(v2, v0);
  mesh.createTriangle(e0, e1, e2);

  Eigen::Vector3d                     center        = Eigen::Vector3d::Constant(0.0);
  Eigen::Vector3d                     halflengths   = Eigen::Vector3d::Constant(1.0);
  FindVoxelContent::BoundaryInclusion includeBounds = FindVoxelContent::INCLUDE_BOUNDARY;
  FindVoxelContent::BoundaryInclusion excludeBounds = FindVoxelContent::EXCLUDE_BOUNDARY;
  query::FindVoxelContent             findIncluded(center, halflengths, includeBounds);
  query::FindVoxelContent             findExcluded(center, halflengths, excludeBounds);

  assertion(testDim >= 0);
  assertion(testDim < dim);

  double sign = positive ? 1.0 : -1.0;

  int thirdDimension = -1;
  for (int dim = 0; dim < 3; dim++) {
    if ((dim != testDim) && (dim != secondDimension)) {
      thirdDimension = dim;
      break;
    }
  }

  // Outside
  coords0[testDim]         = sign * 2.0;
  coords1[testDim]         = sign * 3.0;
  coords2[testDim]         = sign * 2.5;
  coords2[secondDimension] = sign * 0.5;
  v0.setCoords(coords0);
  v1.setCoords(coords1);
  v2.setCoords(coords2);
  mesh.computeState();
  findIncluded(mesh);
  findExcluded(mesh);
  int size = findIncluded.content().triangles().size();
  BOOST_TEST(size == 0);
  size = findExcluded.content().triangles().size();
  BOOST_TEST(size == 0);

  // Outside eps vertex
  coords0[testDim]         = sign * (1.0 + 10.0 * math::NUMERICAL_ZERO_DIFFERENCE);
  coords1[testDim]         = sign * 2.0;
  coords2[testDim]         = sign * 1.5;
  coords2[secondDimension] = sign * 0.5;
  v0.setCoords(coords0);
  v1.setCoords(coords1);
  v2.setCoords(coords2);
  mesh.computeState();
  findIncluded(mesh);
  findExcluded(mesh);
  size = findIncluded.content().triangles().size();
  BOOST_TEST(size == 0);
  size = findExcluded.content().triangles().size();
  BOOST_TEST(size == 0);

  // Outside eps edge
  coords0[testDim]         = sign * (1.0 + 10.0 * math::NUMERICAL_ZERO_DIFFERENCE);
  coords1[testDim]         = sign * 2.0;
  coords2[testDim]         = sign * (1.0 + 10.0 * math::NUMERICAL_ZERO_DIFFERENCE);
  coords2[secondDimension] = sign * 0.5;
  v0.setCoords(coords0);
  v1.setCoords(coords1);
  v2.setCoords(coords2);
  mesh.computeState();
  findIncluded(mesh);
  findExcluded(mesh);
  size = findIncluded.content().triangles().size();
  BOOST_TEST(size == 0);
  size = findExcluded.content().triangles().size();
  BOOST_TEST(size == 0);

  // Touching eps vertex
  coords0[testDim]         = sign * (1.0 + math::NUMERICAL_ZERO_DIFFERENCE);
  coords1[testDim]         = sign * 2.0;
  coords2[testDim]         = sign * 1.5;
  coords2[secondDimension] = sign * 0.5;
  v0.setCoords(coords0);
  v1.setCoords(coords1);
  v2.setCoords(coords2);
  mesh.computeState();
  findIncluded(mesh);
  findExcluded(mesh);
  size = findIncluded.content().triangles().size();
  BOOST_TEST(size == 1);
  size = findExcluded.content().triangles().size();
  BOOST_TEST(size == 0);

  // Touching eps edge
  coords0[testDim]         = sign * (1.0 + math::NUMERICAL_ZERO_DIFFERENCE);
  coords1[testDim]         = sign * 2.0;
  coords2[testDim]         = sign * (1.0 + math::NUMERICAL_ZERO_DIFFERENCE);
  coords2[secondDimension] = sign * 0.5;
  v0.setCoords(coords0);
  v1.setCoords(coords1);
  v2.setCoords(coords2);
  mesh.computeState();
  findIncluded(mesh);
  findExcluded(mesh);
  size = findIncluded.content().triangles().size();
  BOOST_TEST(size == 2);
  size = findExcluded.content().triangles().size();
  BOOST_TEST(size == 0);

  // Touching vertex
  coords0[testDim]         = sign * 1.0;
  coords1[testDim]         = sign * 2.0;
  coords2[testDim]         = sign * 1.5;
  coords2[secondDimension] = sign * 0.5;
  v0.setCoords(coords0);
  v1.setCoords(coords1);
  v2.setCoords(coords2);
  mesh.computeState();
  findIncluded(mesh);
  findExcluded(mesh);
  size = findIncluded.content().triangles().size();
  BOOST_TEST(size == 3);
  size = findExcluded.content().triangles().size();
  BOOST_TEST(size == 0);

  // Touching edge
  coords0[testDim]         = sign * 1.0;
  coords1[testDim]         = sign * 2.0;
  coords2[testDim]         = sign * 1.0;
  coords2[secondDimension] = sign * 0.5;
  v0.setCoords(coords0);
  v1.setCoords(coords1);
  v2.setCoords(coords2);
  mesh.computeState();
  findIncluded(mesh);
  findExcluded(mesh);
  size = findIncluded.content().triangles().size();
  BOOST_TEST(size == 4);
  size = findExcluded.content().triangles().size();
  BOOST_TEST(size == 0);

  // Intersecting eps vertex
  coords0[testDim]         = sign * (1.0 - 10.0 * math::NUMERICAL_ZERO_DIFFERENCE);
  coords1[testDim]         = sign * 2.0;
  coords2[testDim]         = sign * 1.5;
  coords2[secondDimension] = sign * 0.5;
  v0.setCoords(coords0);
  v1.setCoords(coords1);
  v2.setCoords(coords2);
  mesh.computeState();
  findIncluded(mesh);
  findExcluded(mesh);
  size = findIncluded.content().triangles().size();
  BOOST_TEST(size == 5);
  size = findExcluded.content().triangles().size();
  BOOST_TEST(size == 1);

  // Intersecting eps edge
  coords0[testDim]         = sign * (1.0 - 10.0 * math::NUMERICAL_ZERO_DIFFERENCE);
  coords1[testDim]         = sign * 2.0;
  coords2[testDim]         = sign * (1.0 - 10.0 * math::NUMERICAL_ZERO_DIFFERENCE);
  coords2[secondDimension] = sign * 0.5;
  v0.setCoords(coords0);
  v1.setCoords(coords1);
  v2.setCoords(coords2);
  mesh.computeState();
  findIncluded(mesh);
  findExcluded(mesh);
  size = findIncluded.content().triangles().size();
  BOOST_TEST(size == 6);
  size = findExcluded.content().triangles().size();
  BOOST_TEST(size == 2);

  // Intersecting vertex
  coords0[testDim]         = sign * 0.8;
  coords1[testDim]         = sign * 2.0;
  coords2[testDim]         = sign * 1.5;
  coords2[secondDimension] = sign * 0.5;
  v0.setCoords(coords0);
  v1.setCoords(coords1);
  v2.setCoords(coords2);
  mesh.computeState();
  findIncluded(mesh);
  findExcluded(mesh);
  size = findIncluded.content().triangles().size();
  BOOST_TEST(size == 7);
  size = findExcluded.content().triangles().size();
  BOOST_TEST(size == 3);

  // Intersecting edge
  coords0[testDim]         = sign * 0.8;
  coords1[testDim]         = sign * 2.0;
  coords2[testDim]         = sign * 0.8;
  coords2[secondDimension] = sign * 0.5;
  v0.setCoords(coords0);
  v1.setCoords(coords1);
  v2.setCoords(coords2);
  mesh.computeState();
  findIncluded(mesh);
  findExcluded(mesh);
  size = findIncluded.content().triangles().size();
  BOOST_TEST(size == 8);
  size = findExcluded.content().triangles().size();
  BOOST_TEST(size == 4);

  // Contained
  coords0[testDim]         = sign * 0.3;
  coords1[testDim]         = sign * 0.8;
  coords2[testDim]         = sign * 0.5;
  coords2[secondDimension] = sign * 0.5;
  v0.setCoords(coords0);
  v1.setCoords(coords1);
  v2.setCoords(coords2);
  mesh.computeState();
  findIncluded(mesh);
  findExcluded(mesh);
  size = findIncluded.content().triangles().size();
  BOOST_TEST(size == 9);
  size = findExcluded.content().triangles().size();
  BOOST_TEST(size == 5);

  // Contained filling
  coords0[testDim]         = sign * -0.9;
  coords1[testDim]         = sign * 0.9;
  coords2[testDim]         = 0.0;
  coords2[secondDimension] = sign * 0.9;
  v0.setCoords(coords0);
  v1.setCoords(coords1);
  v2.setCoords(coords2);
  mesh.computeState();
  findIncluded(mesh);
  findExcluded(mesh);
  size = findIncluded.content().triangles().size();
  BOOST_TEST(size == 10);
  size = findExcluded.content().triangles().size();
  BOOST_TEST(size == 6);

  // Contained cutting
  coords0[testDim]         = sign * -1.5;
  coords1[testDim]         = sign * 1.5;
  coords2[testDim]         = 0.0;
  coords2[secondDimension] = sign * 1.5;
  v0.setCoords(coords0);
  v1.setCoords(coords1);
  v2.setCoords(coords2);
  mesh.computeState();
  findIncluded(mesh);
  findExcluded(mesh);
  size = findIncluded.content().triangles().size();
  BOOST_TEST(size == 11);
  size = findExcluded.content().triangles().size();
  BOOST_TEST(size == 7);

  // Contained cutting wide1
  coords0[testDim]         = sign * -10.0;
  coords1[testDim]         = sign * 10.0;
  coords2[testDim]         = 0.0;
  coords2[secondDimension] = sign * 10.0;
  v0.setCoords(coords0);
  v1.setCoords(coords1);
  v2.setCoords(coords2);
  mesh.computeState();
  findIncluded(mesh);
  findExcluded(mesh);
  size = findIncluded.content().triangles().size();
  BOOST_TEST(size == 12);
  size = findExcluded.content().triangles().size();
  BOOST_TEST(size == 8);

  // Contained cutting wide2
  coords0[testDim]         = sign * -10000.0;
  coords1[testDim]         = sign * 10000.0;
  coords2[testDim]         = 0.0;
  coords2[secondDimension] = sign * 10000.0;
  v0.setCoords(coords0);
  v1.setCoords(coords1);
  v2.setCoords(coords2);
  mesh.computeState();
  findIncluded(mesh);
  findExcluded(mesh);
  size = findIncluded.content().triangles().size();
  BOOST_TEST(size == 13);
  size = findExcluded.content().triangles().size();
  BOOST_TEST(size == 9);

  // Touching contained fully
  coords0[testDim]         = sign * 0.3;
  coords0[secondDimension] = sign * 1.0;
  coords1[testDim]         = sign * 0.8;
  coords1[secondDimension] = sign * 1.0;
  coords2[testDim]         = 0.5;
  coords2[secondDimension] = sign * 1.0;
  coords2[thirdDimension]  = sign * 0.2;
  v0.setCoords(coords0);
  v1.setCoords(coords1);
  v2.setCoords(coords2);
  mesh.computeState();
  findIncluded(mesh);
  findExcluded(mesh);
  size = findIncluded.content().triangles().size();
  BOOST_TEST(size == 14);
  size = findExcluded.content().triangles().size();
  BOOST_TEST(size == 9);

  // Touching fully
  coords0[testDim]         = sign * -1.5;
  coords0[secondDimension] = sign * 1.0;
  coords1[testDim]         = sign * 1.5;
  coords1[secondDimension] = sign * 1.0;
  coords2[testDim]         = 0.0;
  coords2[secondDimension] = sign * 1.0;
  coords2[thirdDimension]  = sign * 1.5;
  v0.setCoords(coords0);
  v1.setCoords(coords1);
  v2.setCoords(coords2);
  mesh.computeState();
  findIncluded(mesh);
  //  if((sign = -1) && (testDim == 1) && (secondDimension == 2) && (thirdDimension == 0)){
  //    INFO("------------------------------ sign = " << sign << ", testDim = " << testDim
  //                   << ", secondDimension = " << secondDimension << ", thirdDimension = " << thirdDimension);
  //  }
  findExcluded(mesh);
  size = findIncluded.content().triangles().size();
  BOOST_TEST(size == 15);
  size = findExcluded.content().triangles().size();
  //if((sign = -1) && (testDim == 1) && (secondDimension == 2) && (thirdDimension == 0)){
  //    INFO("############################## triangles = " << size);
  //}
  //  if(size != 9){
  //    ERROR("Aus die Mausss");
  //  }
  BOOST_TEST(size == 9);

  // Touching fully wide
  coords0[testDim]         = sign * -10.0;
  coords0[secondDimension] = sign * 1.0;
  coords1[testDim]         = sign * 10.0;
  coords1[secondDimension] = sign * 1.0;
  coords2[testDim]         = 0.0;
  coords2[secondDimension] = sign * 1.0;
  coords2[thirdDimension]  = sign * 10.0;
  v0.setCoords(coords0);
  v1.setCoords(coords1);
  v2.setCoords(coords2);
  mesh.computeState();
  findIncluded(mesh);
  findExcluded(mesh);
  size = findIncluded.content().triangles().size();
  BOOST_TEST(size == 16);
  size = findExcluded.content().triangles().size();
  BOOST_TEST(size == 9);

  // Touching fully wide2
  coords0[testDim]         = sign * -1000.0;
  coords0[secondDimension] = sign * 1.0;
  coords1[testDim]         = sign * 1000.0;
  coords1[secondDimension] = sign * 1.0;
  coords2[testDim]         = 0.0;
  coords2[secondDimension] = sign * 1000.0;
  coords2[thirdDimension]  = sign * 1000.0;
  v0.setCoords(coords0);
  v1.setCoords(coords1);
  v2.setCoords(coords2);
  mesh.computeState();
  findIncluded(mesh);
  findExcluded(mesh);
  size = findIncluded.content().triangles().size();
  BOOST_TEST(size == 17);
  size = findExcluded.content().triangles().size();
  BOOST_TEST(size == 9);
}

BOOST_AUTO_TEST_CASE(Vertices)
{
  for (int dim = 2; dim <= 3; dim++) {
    for (int testDim = 0; testDim < dim; testDim++) {
      bool            positiveDirection = true;
      bool            negativeDirection = false;
      Eigen::VectorXd offset            = Eigen::VectorXd::Zero(dim);
      performTestVertices(testDim, positiveDirection, offset);
      performTestVertices(testDim, negativeDirection, offset);
      offset = Eigen::VectorXd::Constant(dim, -1.0 + 10.0 * math::NUMERICAL_ZERO_DIFFERENCE);
      performTestVertices(testDim, positiveDirection, offset);
      performTestVertices(testDim, negativeDirection, offset);
      offset = Eigen::VectorXd::Constant(dim, 1.0 - 10.0 * math::NUMERICAL_ZERO_DIFFERENCE);
      performTestVertices(testDim, positiveDirection, offset);
      performTestVertices(testDim, negativeDirection, offset);
    }
  }
}

BOOST_AUTO_TEST_CASE(Edges)
{
  for (int dim = 2; dim <= 3; dim++) {
    for (int testDim = 0; testDim < dim; testDim++) {
      bool            positiveDirection = true;
      bool            negativeDirection = false;
      Eigen::VectorXd offset            = Eigen::VectorXd::Zero(dim);
      performTestEdges(testDim, positiveDirection, offset);
      performTestEdges(testDim, negativeDirection, offset);
      offset = Eigen::VectorXd::Constant(dim, -1.0 + 10.0 * math::NUMERICAL_ZERO_DIFFERENCE);
      performTestEdges(testDim, positiveDirection, offset);
      performTestEdges(testDim, negativeDirection, offset);
      offset = Eigen::VectorXd::Constant(dim, 1.0 - 10.0 * math::NUMERICAL_ZERO_DIFFERENCE);
      performTestEdges(testDim, positiveDirection, offset);
      performTestEdges(testDim, negativeDirection, offset);
    }
  }
}

BOOST_AUTO_TEST_CASE(testZeroVoxel)
{
  int           dim = 2;
  mesh::Mesh    mesh("Mesh", dim, false);
  mesh::Vertex &v1 = mesh.createVertex(Eigen::Vector2d(2.0, 0.0));
  mesh::Vertex &v2 = mesh.createVertex(Eigen::Vector2d(2.0, 1.0));
  mesh.createEdge(v1, v2);
  Eigen::Vector2d center(1.0, 1.0);
  Eigen::Vector2d halflengths(0.0, 0.0);

  query::FindVoxelContent find(
      center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
  find(mesh);
  size_t numberContained = find.content().size();
  BOOST_TEST(numberContained == 0);

  mesh.createVertex(Eigen::Vector2d(1.0, 1.1));
  mesh.createVertex(Eigen::Vector2d(1.0, 1.0));
  query::FindVoxelContent find2(center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
  find2(mesh);
  numberContained = find2.content().size();
  BOOST_TEST(numberContained == 1);
}

BOOST_AUTO_TEST_CASE(Triangles)
{
  for (int testDim = 0; testDim < 3; testDim++) {
    bool positiveDirection = true;
    bool negativeDirection = false;
    for (int secondDim = 0; secondDim < 3; secondDim++) {
      if (secondDim == testDim) {
        continue;
      }
      performTestTriangles(testDim, secondDim, positiveDirection);
      performTestTriangles(testDim, secondDim, negativeDirection);
    }
  }
}

BOOST_AUTO_TEST_CASE(CompletelyInsideTriangles)
{
  int             dim              = 3;
  Eigen::Vector3d voxelCenter      = Eigen::Vector3d::Zero();
  Eigen::Vector3d voxelHalflengths = Eigen::Vector3d::Constant(1);

  // Test1
  mesh::Mesh    mesh("Mesh", dim, true);
  mesh::Vertex *v0 = &mesh.createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  mesh::Vertex *v1 = &mesh.createVertex(Eigen::Vector3d(0.5, 0.5, 0.5));
  mesh::Vertex *v2 = &mesh.createVertex(Eigen::Vector3d(0.5, 0.0, 0.5));
  mesh::Edge *  e0 = &mesh.createEdge(*v0, *v1);
  mesh::Edge *  e1 = &mesh.createEdge(*v1, *v2);
  mesh::Edge *  e2 = &mesh.createEdge(*v2, *v0);
  mesh.createTriangle(*e0, *e1, *e2);
  mesh.computeState();

  query::FindVoxelContent find(voxelCenter, voxelHalflengths, FindVoxelContent::INCLUDE_BOUNDARY);
  find(mesh);
  size_t count = find.content().triangles().size();
  BOOST_TEST(count == 1);

  // Test 2
  v0 = &mesh.createVertex(Eigen::Vector3d(0.9, 0.0, 0.0));
  v1 = &mesh.createVertex(Eigen::Vector3d(0.9, 0.9, 0.0));
  v2 = &mesh.createVertex(Eigen::Vector3d(0.9, 0.9, 0.9));
  e0 = &mesh.createEdge(*v0, *v1);
  e1 = &mesh.createEdge(*v1, *v2);
  e2 = &mesh.createEdge(*v2, *v0);
  mesh.createTriangle(*e0, *e1, *e2);
  mesh.computeState();

  find.clear();
  find(mesh);
  count = find.content().triangles().size();
  BOOST_TEST(count == 2);

  // Test 3
  v0 = &mesh.createVertex(Eigen::Vector3d(-0.9, 0.0, 0.0));
  v1 = &mesh.createVertex(Eigen::Vector3d(-0.9, -0.9, 0.0));
  v2 = &mesh.createVertex(Eigen::Vector3d(-0.9, -0.9, -0.9));
  e0 = &mesh.createEdge(*v0, *v1);
  e1 = &mesh.createEdge(*v1, *v2);
  e2 = &mesh.createEdge(*v2, *v0);
  mesh.createTriangle(*e0, *e1, *e2);
  mesh.computeState();

  find.clear();
  find(mesh);
  count = find.content().triangles().size();
  BOOST_TEST(count == 3);
}

BOOST_AUTO_TEST_CASE(CompletelyOutsideTriangles)
{
  int             dim              = 3;
  Eigen::Vector3d voxelCenter      = Eigen::Vector3d::Zero();
  Eigen::Vector3d voxelHalflengths = Eigen::Vector3d::Constant(1.0);

  // Test1
  mesh::Mesh    mesh("Mesh", dim, true);
  mesh::Vertex *v0 = &mesh.createVertex(Eigen::Vector3d(2.0, 2.0, 2.0));
  mesh::Vertex *v1 = &mesh.createVertex(Eigen::Vector3d(2.5, 2.5, 2.5));
  mesh::Vertex *v2 = &mesh.createVertex(Eigen::Vector3d(2.5, 2.0, 2.5));
  mesh::Edge *  e0 = &mesh.createEdge(*v0, *v1);
  mesh::Edge *  e1 = &mesh.createEdge(*v1, *v2);
  mesh::Edge *  e2 = &mesh.createEdge(*v2, *v0);
  mesh.createTriangle(*e0, *e1, *e2);
  mesh.computeState();

  query::FindVoxelContent find(voxelCenter, voxelHalflengths, FindVoxelContent::INCLUDE_BOUNDARY);
  find(mesh);
  size_t count = find.content().triangles().size();
  BOOST_TEST(count == 0);

  // Test 2
  v0 = &mesh.createVertex(Eigen::Vector3d(1.1, 0.0, 0.0));
  v1 = &mesh.createVertex(Eigen::Vector3d(1.1, 1.1, 0.0));
  v2 = &mesh.createVertex(Eigen::Vector3d(1.1, 1.1, 1.1));
  e0 = &mesh.createEdge(*v0, *v1);
  e1 = &mesh.createEdge(*v1, *v2);
  e2 = &mesh.createEdge(*v2, *v0);
  mesh.createTriangle(*e0, *e1, *e2);
  mesh.computeState();

  find.clear();
  find(mesh);
  count = find.content().triangles().size();
  BOOST_TEST(count == 0);

  // Test 3
  v0 = &mesh.createVertex(Eigen::Vector3d(-1.1, 0.0, 0.0));
  v1 = &mesh.createVertex(Eigen::Vector3d(-1.1, -1.1, 0.0));
  v2 = &mesh.createVertex(Eigen::Vector3d(-1.1, -1.1, -1.1));
  e0 = &mesh.createEdge(*v0, *v1);
  e1 = &mesh.createEdge(*v1, *v2);
  e2 = &mesh.createEdge(*v2, *v0);
  mesh.createTriangle(*e0, *e1, *e2);
  mesh.computeState();

  find.clear();
  find(mesh);
  count = find.content().triangles().size();
  BOOST_TEST(count == 0);
}

BOOST_AUTO_TEST_CASE(IntersectingTriangles)
{
  int             dim              = 3;
  Eigen::Vector3d voxelCenter      = Eigen::Vector3d::Zero();
  Eigen::Vector3d voxelHalflengths = Eigen::Vector3d::Constant(1.0);

  // Test1
  mesh::Mesh    mesh("Mesh", dim, true);
  mesh::Vertex *v0 = &mesh.createVertex(Eigen::Vector3d(0.5, 0.0, 0.0));
  mesh::Vertex *v1 = &mesh.createVertex(Eigen::Vector3d(2.0, 1.0, 1.0));
  mesh::Vertex *v2 = &mesh.createVertex(Eigen::Vector3d(2.0, -1.0, -1.0));
  mesh::Edge *  e0 = &mesh.createEdge(*v0, *v1);
  mesh::Edge *  e1 = &mesh.createEdge(*v1, *v2);
  mesh::Edge *  e2 = &mesh.createEdge(*v2, *v0);
  mesh.createTriangle(*e0, *e1, *e2);
  mesh.computeState();

  query::FindVoxelContent find(voxelCenter, voxelHalflengths, FindVoxelContent::INCLUDE_BOUNDARY);
  find(mesh);
  size_t count = find.content().triangles().size();
  BOOST_TEST(count == 1);

  // Test 2
  v0 = &mesh.createVertex(Eigen::Vector3d(-0.5, 0.0, 0.0));
  v1 = &mesh.createVertex(Eigen::Vector3d(-2.0, 1.0, 1.0));
  v2 = &mesh.createVertex(Eigen::Vector3d(-2.0, -1.0, -1.0));
  e0 = &mesh.createEdge(*v0, *v1);
  e1 = &mesh.createEdge(*v1, *v2);
  e2 = &mesh.createEdge(*v2, *v0);
  mesh.createTriangle(*e0, *e1, *e2);
  mesh.computeState();

  find.clear();
  find(mesh);
  count = find.content().triangles().size();
  BOOST_TEST(count == 2);

  // Test 3
  v0 = &mesh.createVertex(Eigen::Vector3d(0.0, 0.5, 0.0));
  v1 = &mesh.createVertex(Eigen::Vector3d(1.0, 2.0, 1.0));
  v2 = &mesh.createVertex(Eigen::Vector3d(-1.0, 2.0, -1.0));
  e0 = &mesh.createEdge(*v0, *v1);
  e1 = &mesh.createEdge(*v1, *v2);
  e2 = &mesh.createEdge(*v2, *v0);
  mesh.createTriangle(*e0, *e1, *e2);
  mesh.computeState();

  find.clear();
  find(mesh);
  count = find.content().triangles().size();
  BOOST_TEST(count == 3);

  // Test 4
  v0 = &mesh.createVertex(Eigen::Vector3d(0.0, -0.5, 0.0));
  v1 = &mesh.createVertex(Eigen::Vector3d(1.0, -2.0, 1.0));
  v2 = &mesh.createVertex(Eigen::Vector3d(-1.0, -2.0, -1.0));
  e0 = &mesh.createEdge(*v0, *v1);
  e1 = &mesh.createEdge(*v1, *v2);
  e2 = &mesh.createEdge(*v2, *v0);
  mesh.createTriangle(*e0, *e1, *e2);
  mesh.computeState();

  find.clear();
  find(mesh);
  count = find.content().triangles().size();
  BOOST_TEST(count == 4);

  // Test 5
  v0 = &mesh.createVertex(Eigen::Vector3d(0.0, 0.0, 0.5));
  v1 = &mesh.createVertex(Eigen::Vector3d(1.0, 1.0, 2.0));
  v2 = &mesh.createVertex(Eigen::Vector3d(-1.0, -1.0, 2.0));
  e0 = &mesh.createEdge(*v0, *v1);
  e1 = &mesh.createEdge(*v1, *v2);
  e2 = &mesh.createEdge(*v2, *v0);
  mesh.createTriangle(*e0, *e1, *e2);
  mesh.computeState();

  find.clear();
  find(mesh);
  count = find.content().triangles().size();
  BOOST_TEST(count == 5);

  // Test 6
  v0 = &mesh.createVertex(Eigen::Vector3d(0.0, 0.0, -0.5));
  v1 = &mesh.createVertex(Eigen::Vector3d(1.0, 1.0, -2.0));
  v2 = &mesh.createVertex(Eigen::Vector3d(-1.0, -1.0, -2.0));
  e0 = &mesh.createEdge(*v0, *v1);
  e1 = &mesh.createEdge(*v1, *v2);
  e2 = &mesh.createEdge(*v2, *v0);
  mesh.createTriangle(*e0, *e1, *e2);
  mesh.computeState();

  find.clear();
  find(mesh);
  count = find.content().triangles().size();
  BOOST_TEST(count == 6);

  // Test 6, triangle cuts out corner of voxel, no triangle vertices contained
  v0 = &mesh.createVertex(Eigen::Vector3d(0.5, 1.5, -0.0));
  v1 = &mesh.createVertex(Eigen::Vector3d(0.5, 0.5, -1.5));
  v2 = &mesh.createVertex(Eigen::Vector3d(1.5, 0.5, -0.5));
  e0 = &mesh.createEdge(*v0, *v1);
  e1 = &mesh.createEdge(*v1, *v2);
  e2 = &mesh.createEdge(*v2, *v0);
  mesh.createTriangle(*e0, *e1, *e2);
  mesh.computeState();

  find.clear();
  find(mesh);
  count = find.content().triangles().size();
  BOOST_TEST(count == 7);
}

BOOST_AUTO_TEST_CASE(TouchingTriangles)
{
  int             dim              = 3;
  Eigen::Vector3d voxelCenter      = Eigen::Vector3d::Zero();
  Eigen::Vector3d voxelHalflengths = Eigen::Vector3d::Constant(1.0);

  // Test1, first triangle vertex touches voxel side
  mesh::Mesh    mesh("Mesh", dim, true);
  mesh::Vertex *v0 = &mesh.createVertex(Eigen::Vector3d(1.0, 0.0, 0.0));
  mesh::Vertex *v1 = &mesh.createVertex(Eigen::Vector3d(2.0, 1.0, 1.0));
  mesh::Vertex *v2 = &mesh.createVertex(Eigen::Vector3d(2.0, -1.0, -1.0));
  mesh::Edge *  e0 = &mesh.createEdge(*v0, *v1);
  mesh::Edge *  e1 = &mesh.createEdge(*v1, *v2);
  mesh::Edge *  e2 = &mesh.createEdge(*v2, *v0);
  mesh.createTriangle(*e0, *e1, *e2);
  mesh.computeState();

  query::FindVoxelContent findIncluded(voxelCenter, voxelHalflengths, FindVoxelContent::INCLUDE_BOUNDARY);
  query::FindVoxelContent findExcluded(voxelCenter, voxelHalflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
  findIncluded(mesh);
  size_t count = findIncluded.content().triangles().size();
  BOOST_TEST(count == 1);
  findExcluded(mesh);
  count = findExcluded.content().triangles().size();
  BOOST_TEST(count == 0);

  // Test2, second triangle vertex touches voxel side
  v0 = &mesh.createVertex(Eigen::Vector3d(1.0, -2.0, 1.0));
  v1 = &mesh.createVertex(Eigen::Vector3d(0.0, -1.0, 0.0));
  v2 = &mesh.createVertex(Eigen::Vector3d(1.0, -2.0, -1.0));
  e0 = &mesh.createEdge(*v0, *v1);
  e1 = &mesh.createEdge(*v1, *v2);
  e2 = &mesh.createEdge(*v2, *v0);
  mesh.createTriangle(*e0, *e1, *e2);
  mesh.computeState();

  findIncluded.clear();
  findIncluded(mesh);
  count = findIncluded.content().triangles().size();
  BOOST_TEST(count == 2);
  findExcluded.clear();
  findExcluded(mesh);
  count = findExcluded.content().triangles().size();
  BOOST_TEST(count == 0);

  // Test3, third triangle vertex touches voxel side
  v0 = &mesh.createVertex(Eigen::Vector3d(-1.0, 0.0, 2.0));
  v1 = &mesh.createVertex(Eigen::Vector3d(1.0, 0.0, 2.0));
  v2 = &mesh.createVertex(Eigen::Vector3d(0.0, 0.0, 1.0));
  e0 = &mesh.createEdge(*v0, *v1);
  e1 = &mesh.createEdge(*v1, *v2);
  e2 = &mesh.createEdge(*v2, *v0);
  mesh.createTriangle(*e0, *e1, *e2);
  mesh.computeState();

  findIncluded.clear();
  findIncluded(mesh);
  count = findIncluded.content().triangles().size();
  BOOST_TEST(count == 3);
  findExcluded.clear();
  findExcluded(mesh);
  count = findExcluded.content().triangles().size();
  BOOST_TEST(count == 0);

  // Test4, triangle edge touches voxel side
  v0 = &mesh.createVertex(Eigen::Vector3d(1.0, -1.0, -1.0));
  v1 = &mesh.createVertex(Eigen::Vector3d(1.0, 1.0, 1.0));
  v2 = &mesh.createVertex(Eigen::Vector3d(2.0, 0.0, 0.0));
  e0 = &mesh.createEdge(*v0, *v1);
  e1 = &mesh.createEdge(*v1, *v2);
  e2 = &mesh.createEdge(*v2, *v0);
  mesh.createTriangle(*e0, *e1, *e2);
  mesh.computeState();

  findIncluded.clear();
  findIncluded(mesh);
  count = findIncluded.content().triangles().size();
  BOOST_TEST(count == 4);
  findExcluded.clear();
  findExcluded(mesh);
  count = findExcluded.content().triangles().size();
  BOOST_TEST(count == 0);

  // Test5, triangle edge touches voxel edge in one point
  v0 = &mesh.createVertex(Eigen::Vector3d(2.0, 0.0, 0.0));
  v1 = &mesh.createVertex(Eigen::Vector3d(0.0, 2.0, 0.0));
  v2 = &mesh.createVertex(Eigen::Vector3d(2.0, 2.0, 0.0));
  e0 = &mesh.createEdge(*v0, *v1);
  e1 = &mesh.createEdge(*v1, *v2);
  e2 = &mesh.createEdge(*v2, *v0);
  mesh.createTriangle(*e0, *e1, *e2);
  mesh.computeState();

  findIncluded.clear();
  findIncluded(mesh);
  count = findIncluded.content().triangles().size();
  BOOST_TEST(count == 5);
  findExcluded.clear();
  findExcluded(mesh);
  count = findExcluded.content().triangles().size();
  BOOST_TEST(count == 0);

  // Test6, triangle edge overlays with voxel edge
  v0 = &mesh.createVertex(Eigen::Vector3d(-2.0, -1.0, 1.0));
  v1 = &mesh.createVertex(Eigen::Vector3d(0.0, -1.0, 1.0));
  v2 = &mesh.createVertex(Eigen::Vector3d(-1.0, 2.0, 2.0));
  e0 = &mesh.createEdge(*v0, *v1);
  e1 = &mesh.createEdge(*v1, *v2);
  e2 = &mesh.createEdge(*v2, *v0);
  mesh.createTriangle(*e0, *e1, *e2);
  mesh.computeState();

  findIncluded.clear();
  findIncluded(mesh);
  count = findIncluded.content().triangles().size();
  BOOST_TEST(count == 6);
  findExcluded.clear();
  findExcluded(mesh);
  count = findExcluded.content().triangles().size();
  BOOST_TEST(count == 0);

  // Test6, triangle is contained in voxel side
  v0 = &mesh.createVertex(Eigen::Vector3d(-0.5, -1.0, 0.0));
  v1 = &mesh.createVertex(Eigen::Vector3d(0.5, -1.0, 0.0));
  v2 = &mesh.createVertex(Eigen::Vector3d(0.5, -1.0, 0.5));
  e0 = &mesh.createEdge(*v0, *v1);
  e1 = &mesh.createEdge(*v1, *v2);
  e2 = &mesh.createEdge(*v2, *v0);
  mesh.createTriangle(*e0, *e1, *e2);
  mesh.computeState();

  findIncluded.clear();
  findIncluded(mesh);
  count = findIncluded.content().triangles().size();
  BOOST_TEST(count == 7);
  findExcluded.clear();
  findExcluded(mesh);
  count = findExcluded.content().triangles().size();
  BOOST_TEST(count == 0);

  // Test6, voxel side is contained in triangle
  v0 = &mesh.createVertex(Eigen::Vector3d(-5.0, -5.0, 1.0));
  v1 = &mesh.createVertex(Eigen::Vector3d(5.0, -5.0, 1.0));
  v2 = &mesh.createVertex(Eigen::Vector3d(0.0, 5.0, 1.0));
  e0 = &mesh.createEdge(*v0, *v1);
  e1 = &mesh.createEdge(*v1, *v2);
  e2 = &mesh.createEdge(*v2, *v0);
  mesh.createTriangle(*e0, *e1, *e2);
  mesh.computeState();

  findIncluded.clear();
  findIncluded(mesh);
  count = findIncluded.content().triangles().size();
  BOOST_TEST(count == 8);
  findExcluded.clear();
  findExcluded(mesh);
  count = findExcluded.content().triangles().size();
  BOOST_TEST(count == 0);
}

BOOST_AUTO_TEST_CASE(QueryCube)
{
  using namespace mesh;
  int     dim         = 3;
  bool    flipNormals = false;
  Mesh    mesh("TestMesh", dim, flipNormals);
  Vertex &v000 = mesh.createVertex(Eigen::Vector3d(0.0, 0.0, 0.0));
  Vertex &v001 = mesh.createVertex(Eigen::Vector3d(0.0, 0.0, 1.0));
  Vertex &v010 = mesh.createVertex(Eigen::Vector3d(0.0, 1.0, 0.0));
  Vertex &v011 = mesh.createVertex(Eigen::Vector3d(0.0, 1.0, 1.0));
  Vertex &v100 = mesh.createVertex(Eigen::Vector3d(1.0, 0.0, 0.0));
  Vertex &v101 = mesh.createVertex(Eigen::Vector3d(1.0, 0.0, 1.0));
  Vertex &v110 = mesh.createVertex(Eigen::Vector3d(1.0, 1.0, 0.0));
  Vertex &v111 = mesh.createVertex(Eigen::Vector3d(1.0, 1.0, 1.0));

  Edge &e000to100 = mesh.createEdge(v000, v100);
  Edge &e010to110 = mesh.createEdge(v010, v110);
  Edge &e001to101 = mesh.createEdge(v001, v101);
  Edge &e011to111 = mesh.createEdge(v011, v111);

  Edge &e000to010 = mesh.createEdge(v000, v010);
  Edge &e100to110 = mesh.createEdge(v100, v110);
  Edge &e001to011 = mesh.createEdge(v001, v011);
  Edge &e101to111 = mesh.createEdge(v101, v111);

  Edge &e000to001 = mesh.createEdge(v000, v001);
  Edge &e100to101 = mesh.createEdge(v100, v101);
  Edge &e010to011 = mesh.createEdge(v010, v011);
  Edge &e110to111 = mesh.createEdge(v110, v111);

  Edge &e010to001 = mesh.createEdge(v010, v001);

  Edge &e100to111 = mesh.createEdge(v100, v111);
  Edge &e100to001 = mesh.createEdge(v100, v001);
  Edge &e110to011 = mesh.createEdge(v110, v011);

  Edge &e010to100 = mesh.createEdge(v010, v100);

  Edge &e001to111 = mesh.createEdge(v001, v111);

  mesh.createTriangle(e000to001, e010to001, e000to010); // x = 0
  mesh.createTriangle(e010to011, e010to001, e001to011);
  mesh.createTriangle(e100to101, e100to111, e101to111); // x = 1
  mesh.createTriangle(e100to110, e110to111, e100to111);

  mesh.createTriangle(e000to100, e100to001, e000to001); // y = 0
  mesh.createTriangle(e100to101, e001to101, e100to001);
  mesh.createTriangle(e010to110, e010to011, e110to011); // y = 1
  mesh.createTriangle(e110to111, e110to011, e011to111);

  mesh.createTriangle(e010to100, e000to100, e000to010); // z = 0
  mesh.createTriangle(e010to100, e010to110, e100to110);
  mesh.createTriangle(e001to101, e101to111, e001to111); // z = 1
  mesh.createTriangle(e001to011, e001to111, e011to111);

  mesh.computeState();

  io::ExportVTK exportMesh(true);
  std::string   location = "";
  exportMesh.doExport("FindVoxelContentTest-testQueryCube", location, mesh);

  // Query mesh
  Eigen::Vector3d center      = Eigen::Vector3d::Zero();
  Eigen::Vector3d halflengths = Eigen::Vector3d::Zero();

  // Vertex queries with halflengths = 1.0/3.0
  {
    center << 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0;
    halflengths = Eigen::Vector3d::Constant(1.0 / 3.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 2.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0;
    halflengths = Eigen::Vector3d::Constant(1.0 / 3.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 1.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0;
    halflengths = Eigen::Vector3d::Constant(1.0 / 3.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 2.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0;
    halflengths = Eigen::Vector3d::Constant(1.0 / 3.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 1.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0;
    halflengths = Eigen::Vector3d::Constant(1.0 / 3.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 2.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0;
    halflengths = Eigen::Vector3d::Constant(1.0 / 3.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0;
    halflengths = Eigen::Vector3d::Constant(1.0 / 3.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 2.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0;
    halflengths = Eigen::Vector3d::Constant(1.0 / 3.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }

  // z = 1.0/6.0 plane
  {
    center << 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 0.5, 1.0 / 6.0, 1.0 / 6.0;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 5.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 1.0 / 6.0, 0.5, 1.0 / 6.0;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 0.5, 0.5, 1.0 / 6.0;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 5.0 / 6.0, 0.5, 1.0 / 6.0;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 1.0 / 6.0, 5.0 / 6.0, 1.0 / 6.0;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 0.5, 5.0 / 6.0, 1.0 / 6.0;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 5.0 / 6.0, 5.0 / 6.0, 1.0 / 6.0;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }

  // z = 0.5 plane
  {
    center << 1.0 / 6.0, 1.0 / 6.0, 0.5;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 0.5, 1.0 / 6.0, 0.5;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 5.0 / 6.0, 1.0 / 6.0, 0.5;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 1.0 / 6.0, 0.5, 0.5;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 0.5, 0.5, 0.5;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() == 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 5.0 / 6.0, 0.5, 0.5;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 1.0 / 6.0, 5.0 / 6.0, 0.5;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 0.5, 5.0 / 6.0, 0.5;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 5.0 / 6.0, 5.0 / 6.0, 0.5;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }

  // z = 5.0/6.0 plane
  {
    center << 1.0 / 6.0, 1.0 / 6.0, 5.0 / 6.0;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 0.5, 1.0 / 6.0, 5.0 / 6.0;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 5.0 / 6.0, 1.0 / 6.0, 5.0 / 6.0;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 1.0 / 6.0, 0.5, 5.0 / 6.0;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 0.5, 0.5, 5.0 / 6.0;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 5.0 / 6.0, 0.5, 5.0 / 6.0;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 1.0 / 6.0, 5.0 / 6.0, 5.0 / 6.0;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 0.5, 5.0 / 6.0, 5.0 / 6.0;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 5.0 / 6.0, 5.0 / 6.0, 5.0 / 6.0;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }

  // Query voxel equals mesh outline
  {
    center << 0.5, 0.5, 0.5;
    halflengths = Eigen::Vector3d::Constant(0.5);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().vertices().size() == 8);
    BOOST_TEST(findIncluded.content().edges().size() == 18);
    BOOST_TEST(findIncluded.content().triangles().size() == 12);
    BOOST_TEST(findExcluded.content().size() == 0);
  }

  // Invert mesh normal and perform all queries again, results have to be same
  flipNormals = true;
  mesh.setFlipNormals(flipNormals);
  mesh.computeState();

  // Vertex queries with halflengths = 1.0/3.0
  {
    center << 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0;
    halflengths = Eigen::Vector3d::Constant(1.0 / 3.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 2.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0;
    halflengths = Eigen::Vector3d::Constant(1.0 / 3.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 1.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0;
    halflengths = Eigen::Vector3d::Constant(1.0 / 3.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 2.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0;
    halflengths = Eigen::Vector3d::Constant(1.0 / 3.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 1.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0;
    halflengths = Eigen::Vector3d::Constant(1.0 / 3.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 2.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0;
    halflengths = Eigen::Vector3d::Constant(1.0 / 3.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0;
    halflengths = Eigen::Vector3d::Constant(1.0 / 3.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 2.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0;
    halflengths = Eigen::Vector3d::Constant(1.0 / 3.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }

  // z = 1.0/6.0 plane
  {
    center << 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 0.5, 1.0 / 6.0, 1.0 / 6.0;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 5.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 1.0 / 6.0, 0.5, 1.0 / 6.0;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 0.5, 0.5, 1.0 / 6.0;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 5.0 / 6.0, 0.5, 1.0 / 6.0;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 1.0 / 6.0, 5.0 / 6.0, 1.0 / 6.0;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 0.5, 5.0 / 6.0, 1.0 / 6.0;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 5.0 / 6.0, 5.0 / 6.0, 1.0 / 6.0;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }

  // z = 0.5 plane
  {
    center << 1.0 / 6.0, 1.0 / 6.0, 0.5;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 0.5, 1.0 / 6.0, 0.5;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 5.0 / 6.0, 1.0 / 6.0, 0.5;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 1.0 / 6.0, 0.5, 0.5;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 0.5, 0.5, 0.5;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() == 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 5.0 / 6.0, 0.5, 0.5;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 1.0 / 6.0, 5.0 / 6.0, 0.5;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 0.5, 5.0 / 6.0, 0.5;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 5.0 / 6.0, 5.0 / 6.0, 0.5;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }

  // z = 5.0/6.0 plane
  {
    center << 1.0 / 6.0, 1.0 / 6.0, 5.0 / 6.0;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 0.5, 1.0 / 6.0, 5.0 / 6.0;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 5.0 / 6.0, 1.0 / 6.0, 5.0 / 6.0;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 1.0 / 6.0, 0.5, 5.0 / 6.0;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 0.5, 0.5, 5.0 / 6.0;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 5.0 / 6.0, 0.5, 5.0 / 6.0;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 1.0 / 6.0, 5.0 / 6.0, 5.0 / 6.0;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 0.5, 5.0 / 6.0, 5.0 / 6.0;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }
  {
    center << 5.0 / 6.0, 5.0 / 6.0, 5.0 / 6.0;
    halflengths = Eigen::Vector3d::Constant(1.0 / 6.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().size() == 0);
  }

  // Query voxel equals mesh outline
  {
    center << 0.5, 0.5, 0.5;
    halflengths = Eigen::Vector3d::Constant(0.5);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    BOOST_TEST(findIncluded.content().vertices().size() == 8);
    BOOST_TEST(findIncluded.content().edges().size() == 18);
    BOOST_TEST(findIncluded.content().triangles().size() == 12);
    BOOST_TEST(findExcluded.content().size() == 0);
  }

  // Special query that gave error in real scenario
  {
    center << 1.0 / 9.0, 8.0 / 9.0, -5.5511151231257827e-17;
    halflengths = Eigen::Vector3d::Constant(1.0 / 9.0);
    query::FindVoxelContent findIncluded(
        center, halflengths, FindVoxelContent::INCLUDE_BOUNDARY);
    query::FindVoxelContent findExcluded(
        center, halflengths, FindVoxelContent::EXCLUDE_BOUNDARY);
    findIncluded(mesh);
    findExcluded(mesh);
    ExportVTKVoxelQueries exportQueries;
    exportQueries.addQuery(center, halflengths, 0);
    exportQueries.exportQueries("FindVoxelContentTest-testQueryCube-queries");
    BOOST_TEST(findIncluded.content().size() > 0);
    BOOST_TEST(findExcluded.content().edges().size() == 1);
    BOOST_TEST(findExcluded.content().triangles().size() == 2);
  }
}

BOOST_AUTO_TEST_SUITE_END() // FindVoxelContentTests
BOOST_AUTO_TEST_SUITE_END() // QueryTests
