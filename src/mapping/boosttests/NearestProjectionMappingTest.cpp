#include "testing/Testing.hpp"

#include "mapping/NearestProjectionMapping.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Edge.hpp"

using namespace precice;

BOOST_AUTO_TEST_SUITE(MappingTests)
BOOST_AUTO_TEST_SUITE(NearestProjectionMapping)

BOOST_AUTO_TEST_CASE(testConservativeNonIncremental)                
{
  using namespace mesh;
  int dimensions = 2;

  // Setup geometry to map to
  PtrMesh outMesh ( new Mesh("OutMesh", dimensions, true) );
  PtrData outData = outMesh->createData ( "Data", 1 );
  int outDataID = outData->getID();
  Vertex& v1 = outMesh->createVertex ( Eigen::Vector2d(0.0, 0.0) );
  Vertex& v2 = outMesh->createVertex ( Eigen::Vector2d(1.0, 1.0) );
  outMesh->createEdge ( v1, v2 );
  outMesh->computeState();
  outMesh->allocateDataValues();

  PtrMesh inMesh ( new Mesh("InMesh", dimensions, false) );
  PtrData inData = inMesh->createData ( "Data", 1 );
  int inDataID = inData->getID();

  // Setup mapping with mapping coordinates and geometry used
  mapping::NearestProjectionMapping mapping(mapping::Mapping::CONSERVATIVE, dimensions);
  mapping.setMeshes(inMesh, outMesh);

  // Map value 1.0 from middle of edge to geometry. Expect half of the
  // value to be added to vertex1 and half of it to vertex2.
  Vertex& inv1 = inMesh->createVertex(Eigen::Vector2d(0.5, 0.5));
  // Map value 1.0 from below edge to geometry. Expect vertex1 to get the
  // full data value, i.e. 1.0 and in addition the value from before. In total
  // v1 should have 1.5 * dataValue then.
  Vertex& inv2 = inMesh->createVertex(Eigen::Vector2d(-0.5, -0.5));
  // Do the same thing from above, expect vertex2 to get the full value now.
  Vertex& inv3 = inMesh->createVertex(Eigen::Vector2d(1.5, 1.5));
  inMesh->allocateDataValues();

  double value = 1.0;
  //assign(inData->values()) = value;
  inData->values() = Eigen::VectorXd::Constant(inData->values().size(), value);
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  Eigen::VectorXd& values = outData->values();
  BOOST_TEST(values(0) == value * 1.5);
  BOOST_TEST(values(1) == value * 1.5);

  // Change in-vertex coordinates and recompute mapping
  inv1.setCoords (Eigen::Vector2d(-1.0, -1.0));
  inv2.setCoords (Eigen::Vector2d(-1.0, -1.0));
  inv3.setCoords (Eigen::Vector2d(1.0, 1.0));
  //assign(values) = 0.0;
  values = Eigen::VectorXd::Constant(values.size(), 0.0);

  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(values(0) == value * 2.0);
  BOOST_TEST(values(1) == value * 1.0);

  // reset output value and remap
  //assign(values) = 0.0;
  values = Eigen::VectorXd::Constant(values.size(), 0.0);

  mapping.map(inDataID, outDataID);
  BOOST_TEST(values(0) == value * 2.0);
  BOOST_TEST(values(1) == value * 1.0);
}

BOOST_AUTO_TEST_CASE(ConsistentNonIncremental2D)
{
  using namespace mesh;
  int dimensions = 2;

  // Create mesh to map from
  PtrMesh inMesh ( new Mesh("InMesh", dimensions, false) );
  PtrData inData = inMesh->createData ( "InData", 1 );
  int inDataID = inData->getID ();
  Vertex& v1 = inMesh->createVertex ( Eigen::Vector2d(0.0, 0.0) );
  Vertex& v2 = inMesh->createVertex ( Eigen::Vector2d(1.0, 1.0) );
  inMesh->createEdge ( v1, v2 );
  inMesh->computeState();
  inMesh->allocateDataValues();
  double valueVertex1 = 1.0;
  double valueVertex2 = 2.0;
  Eigen::VectorXd& values = inData->values();
  values(0) = valueVertex1;
  values(1) = valueVertex2;

  // Create mesh to map to
  PtrMesh outMesh ( new Mesh("OutMesh", dimensions, false) );
  PtrData outData = outMesh->createData ( "OutData", 1 );
  int outDataID = outData->getID();

  // Setup mapping with mapping coordinates and geometry used
  mapping::NearestProjectionMapping mapping(mapping::Mapping::CONSISTENT, dimensions);
  mapping.setMeshes ( inMesh, outMesh );
  BOOST_TEST ( mapping.hasComputedMapping() == false );

  Vertex& outv0 = outMesh->createVertex ( Eigen::Vector2d(0.5, 0.5) );
  Vertex& outv1 = outMesh->createVertex ( Eigen::Vector2d(-0.5, -0.5) );
  Vertex& outv2 = outMesh->createVertex ( Eigen::Vector2d(1.5, 1.5) );
  outMesh->allocateDataValues();

  // Compute and perform mapping
  mapping.computeMapping();
  mapping.map ( inDataID, outDataID );

  // Validate results
  BOOST_TEST ( mapping.hasComputedMapping() == true );
  BOOST_TEST ( outData->values()[0] == (valueVertex1 + valueVertex2) * 0.5 );
  BOOST_TEST ( outData->values()[1] == valueVertex1 );
  BOOST_TEST ( outData->values()[2] == valueVertex2 );

  // Redo mapping, results should be
  //assign(outData->values()) = 0.0;
  outData->values() = Eigen::VectorXd::Constant(outData->values().size(), 0.0);

  mapping.map ( inDataID, outDataID );
  BOOST_TEST ( outData->values()[0] == (valueVertex1 + valueVertex2) * 0.5 );
  BOOST_TEST ( outData->values()[1] == valueVertex1 );
  BOOST_TEST ( outData->values()[2] == valueVertex2 );

  // Change vertex coordinates and redo mapping
  outv0.setCoords ( Eigen::Vector2d(-0.5, -0.5) );
  outv1.setCoords ( Eigen::Vector2d(1.5, 1.5) );
  outv2.setCoords ( Eigen::Vector2d(0.5, 0.5) );
  //assign(outData->values()) = 0.0;
  outData->values() = Eigen::VectorXd::Constant(outData->values().size(), 0.0);

  mapping.computeMapping();
  mapping.map ( inDataID, outDataID );
  BOOST_TEST ( outData->values()[0] == valueVertex1 );
  BOOST_TEST ( outData->values()[1] == valueVertex2 );
  BOOST_TEST ( outData->values()[2] == (valueVertex1 + valueVertex2) * 0.5 );

  // Reset output data to zero and redo the mapping
  //assign(outData->values()) = 0.0;
  outData->values() = Eigen::VectorXd::Constant(outData->values().size(), 0.0);

  mapping.map ( inDataID, outDataID );
  BOOST_TEST ( outData->values()[0] == valueVertex1 );
  BOOST_TEST ( outData->values()[1] == valueVertex2 );
  BOOST_TEST ( outData->values()[2] == (valueVertex1 + valueVertex2) * 0.5 );
}


BOOST_AUTO_TEST_CASE(ConsistentNonIncrementalPseudo3D)
{
  using namespace mesh;
  int dimensions = 3;

  // Create mesh to map from
  PtrMesh inMesh ( new Mesh("InMesh", dimensions, false) );
  PtrData inData = inMesh->createData ( "InData", 1 );
  int inDataID = inData->getID ();
  Vertex& v1 = inMesh->createVertex ( Eigen::Vector3d(0.0, 0.0, 0.0) );
  Vertex& v2 = inMesh->createVertex ( Eigen::Vector3d(1.0, 1.0, 0.0) );
  inMesh->createEdge ( v1, v2 );
  inMesh->computeState();
  inMesh->allocateDataValues();
  double valueVertex1 = 1.0;
  double valueVertex2 = 2.0;
  Eigen::VectorXd& values = inData->values();
  values(0) = valueVertex1;
  values(1) = valueVertex2;

  // Create mesh to map to
  PtrMesh outMesh ( new Mesh("OutMesh", dimensions, false) );
  PtrData outData = outMesh->createData ( "OutData", 1 );
  int outDataID = outData->getID();

  // Setup mapping with mapping coordinates and geometry used
  mapping::NearestProjectionMapping mapping(mapping::Mapping::CONSISTENT, dimensions);
  mapping.setMeshes ( inMesh, outMesh );
  BOOST_TEST ( mapping.hasComputedMapping() == false );

  Vertex& outv0 = outMesh->createVertex ( Eigen::Vector3d(0.5, 0.5, 0.0) );
  Vertex& outv1 = outMesh->createVertex ( Eigen::Vector3d(-0.5, -0.5, 0.0) );
  Vertex& outv2 = outMesh->createVertex ( Eigen::Vector3d(1.5, 1.5, 0.0) );
  outMesh->allocateDataValues();

  // Compute and perform mapping
  mapping.computeMapping();
  mapping.map ( inDataID, outDataID );

  // Validate results
  BOOST_TEST ( mapping.hasComputedMapping() == true );
  BOOST_TEST ( outData->values()[0] == (valueVertex1 + valueVertex2) * 0.5 );
  BOOST_TEST ( outData->values()[1] == valueVertex1 );
  BOOST_TEST ( outData->values()[2] == valueVertex2 );

  // Redo mapping, results should be
  //assign(outData->values()) = 0.0;
  outData->values() = Eigen::VectorXd::Constant(outData->values().size(), 0.0);

  mapping.map ( inDataID, outDataID );
  BOOST_TEST ( outData->values()[0] == (valueVertex1 + valueVertex2) * 0.5 );
  BOOST_TEST ( outData->values()[1] == valueVertex1 );
  BOOST_TEST ( outData->values()[2] == valueVertex2 );

  // Change vertex coordinates and redo mapping
  outv0.setCoords ( Eigen::Vector3d(-0.5, -0.5, 0.0) );
  outv1.setCoords ( Eigen::Vector3d(1.5, 1.5, 0.0) );
  outv2.setCoords ( Eigen::Vector3d(0.5, 0.5, 0.0) );
  //assign(outData->values()) = 0.0;
  outData->values() = Eigen::VectorXd::Constant(outData->values().size(), 0.0);

  mapping.computeMapping();
  mapping.map ( inDataID, outDataID );
  BOOST_TEST ( outData->values()[0] == valueVertex1 );
  BOOST_TEST ( outData->values()[1] == valueVertex2 );
  BOOST_TEST ( outData->values()[2] == (valueVertex1 + valueVertex2) * 0.5 );

  // Reset output data to zero and redo the mapping
  //assign(outData->values()) = 0.0;
  outData->values() = Eigen::VectorXd::Constant(outData->values().size(), 0.0);

  mapping.map ( inDataID, outDataID );
  BOOST_TEST ( outData->values()[0] == valueVertex1 );
  BOOST_TEST ( outData->values()[1] == valueVertex2 );
  BOOST_TEST ( outData->values()[2] == (valueVertex1 + valueVertex2) * 0.5 );
}


BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
