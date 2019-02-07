#include "testing/Testing.hpp"

#include "mapping/AxialGeoMultiscaleMapping.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"

using namespace precice;

BOOST_AUTO_TEST_SUITE(MappingTests)
BOOST_AUTO_TEST_SUITE(AxialGeoMultiscaleMapping)

BOOST_AUTO_TEST_CASE(testConsistentSpread)
{
  using namespace mesh;
  int dimensions = 3;

  // Setup geometry to map from
  PtrMesh inMesh ( new Mesh("InMesh", dimensions, false) );
  PtrData inData = inMesh->createData ( "Data", dimensions );
  int inDataID = inData->getID();
  inMesh->createVertex ( Eigen::Vector3d(1.0, 4.0, 1.0) );
  inMesh->allocateDataValues();
  double radius = 2.0;


  // Setup geometry to map to
  PtrMesh outMesh ( new Mesh("OutMesh", dimensions, true) );
  PtrData outData = outMesh->createData ( "Data", dimensions );
  int outDataID = outData->getID();
  outMesh->createVertex ( Eigen::Vector3d(1.0, 4.0, 1.0) ); // center of parabula
  outMesh->createVertex ( Eigen::Vector3d(1.0, 2.0, 1.0) ); // distance of r, so should be 0
  outMesh->createVertex ( Eigen::Vector3d(1.0, 4.0, 0.0) ); // distance of 1
  outMesh->computeState();
  outMesh->allocateDataValues();

  // Setup mapping with mapping coordinates and geometry used
  mapping::AxialGeoMultiscaleMapping mapping(mapping::Mapping::CONSISTENT, dimensions,
                                             mapping::AxialGeoMultiscaleMapping::SPREAD, radius);
  mapping.setMeshes(inMesh, outMesh);

  Eigen::Vector3d inValues = Eigen::Vector3d(1.0, 2.0, 3.0);
  inData->values() = inValues;
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  Eigen::VectorXd& outValues = outData->values();
  BOOST_TEST(outValues.size() = dimensions * 3);
  BOOST_TEST(outValues[0] == inValues[0] * 2.0);
  BOOST_TEST(outValues[1] == inValues[1] * 2.0);
  BOOST_TEST(outValues[2] == inValues[2] * 2.0);
  BOOST_TEST(outValues[3] == 0.0);
  BOOST_TEST(outValues[4] == 0.0);
  BOOST_TEST(outValues[5] == 0.0);
  BOOST_TEST(outValues[6] == inValues[0] * 3.0 / 2.0);
  BOOST_TEST(outValues[7] == inValues[1] * 3.0 / 2.0);
  BOOST_TEST(outValues[8] == inValues[2] * 3.0 / 2.0);
}


BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
