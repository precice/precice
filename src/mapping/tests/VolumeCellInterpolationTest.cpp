#include <Eigen/Core>
#include <algorithm>
#include <memory>
#include "logging/LogMacros.hpp"
#include "mapping/Mapping.hpp"
#include "mapping/NearestNeighborMapping.hpp"
#include "math/constants.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Utils.hpp"
#include "mesh/Vertex.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

using namespace precice;
using namespace precice::mesh;

BOOST_AUTO_TEST_SUITE(MappingTests)
BOOST_AUTO_TEST_SUITE(VolumeCellInterpolation)

BOOST_AUTO_TEST_CASE(ConsistentNonIncremental)
{
  PRECICE_TEST(1_rank);
  int dimensions = 2;
  using testing::equals;



  PtrMesh inMesh(new Mesh("InMesh", dimensions, testing::nextMeshID()));
  PtrData inDataScalar   = inMesh->createData("InDataScalar", 1, 0_dataID);
  int inDataScalarID = inDataScalar->getID();

  Vertex& inVertex0 = inMesh->createVertex(Eigen::Vector2d(0.0, 0.0));
  Vertex& inVertex1 = inMesh->createVertex(Eigen::Vector2d(1.0, 0.0));
  Vertex& inVertex2 = inMesh->createVertex(Eigen::Vector2d(0.0, 1.0));

  Edge& inEdge0 = inMesh->createEdge(inVertex0, inVertex1);
    Edge& inEdge1 = inMesh->createEdge(inVertex1, inVertex2);
  Edge& inEdge2 = inMesh->createEdge(inVertex2, inVertex0);

  Triangle& inTriangle = inMesh->createTriangle(inEdge0, inEdge1, inEdge2);
  BOOST_CHECK(!inMesh->edges().empty());
  BOOST_CHECK(!inMesh->triangles().empty());
  
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
