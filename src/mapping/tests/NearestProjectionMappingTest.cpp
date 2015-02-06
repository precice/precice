// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "NearestProjectionMappingTest.hpp"
#include "mapping/NearestProjectionMapping.hpp"
#include "utils/Parallel.hpp"
#include "utils/Globals.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Edge.hpp"

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::mapping::tests::NearestProjectionMappingTest)

namespace precice {
namespace mapping {
namespace tests {

tarch::logging::Log NearestProjectionMappingTest::
  _log ("precice::mapping::tests::NearestProjectionMappingTest");

NearestProjectionMappingTest:: NearestProjectionMappingTest()
:
  tarch::tests::TestCase ("precice::mapping::tests::NearestProjectionMappingTest")
{}

void NearestProjectionMappingTest:: run()
{
  PRECICE_MASTER_ONLY {
    testMethod(testConservativeNonIncremental);
    testMethod(testConsistentNonIncremental);
  }
}

void NearestProjectionMappingTest:: testConservativeNonIncremental()
{
  preciceTrace ( "testConservativeNonIncremental()" );
  using namespace mesh;
  using utils::Vector2D;
  int dimensions = 2;

  // Setup geometry to map to
  PtrMesh outMesh ( new Mesh("OutMesh", dimensions, true) );
  PtrData outData = outMesh->createData ( "Data", 1 );
  int outDataID = outData->getID();
  Vertex& v1 = outMesh->createVertex ( Vector2D(0.0, 0.0) );
  Vertex& v2 = outMesh->createVertex ( Vector2D(1.0, 1.0) );
  outMesh->createEdge ( v1, v2 );
  outMesh->computeState();
  outMesh->allocateDataValues();

  PtrMesh inMesh ( new Mesh("InMesh", dimensions, false) );
  PtrData inData = inMesh->createData ( "Data", 1 );
  int inDataID = inData->getID();

  // Setup mapping with mapping coordinates and geometry used
  NearestProjectionMapping mapping(Mapping::CONSERVATIVE, dimensions);
  mapping.setMeshes(inMesh, outMesh);

  // Map value 1.0 from middle of edge to geometry. Expect half of the
  // value to be added to vertex1 and half of it to vertex2.
  Vertex& inv1 = inMesh->createVertex(Vector2D(0.5, 0.5));
  // Map value 1.0 from below edge to geometry. Expect vertex1 to get the
  // full data value, i.e. 1.0 and in addition the value from before. In total
  // v1 should have 1.5 * dataValue then.
  Vertex& inv2 = inMesh->createVertex(Vector2D(-0.5, -0.5));
  // Do the same thing from above, expect vertex2 to get the full value now.
  Vertex& inv3 = inMesh->createVertex(Vector2D(1.5, 1.5));
  inMesh->allocateDataValues();

  double value = 1.0;
  assign(inData->values()) = value;
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  utils::DynVector& values = outData->values();
  validateNumericalEquals(values[0], value * 1.5);
  validateNumericalEquals(values[1], value * 1.5);

  // Change in-vertex coordinates and recompute mapping
  inv1.setCoords (Vector2D(-1.0, -1.0));
  inv2.setCoords (Vector2D(-1.0, -1.0));
  inv3.setCoords (Vector2D(1.0, 1.0));
  assign(values) = 0.0;
  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  validateNumericalEquals(values[0], value * 2.0);
  validateNumericalEquals(values[1], value * 1.0);

  // reset output value and remap
  assign(values) = 0.0;
  mapping.map(inDataID, outDataID);
  validateNumericalEquals(values[0], value * 2.0);
  validateNumericalEquals(values[1], value * 1.0);
}

void NearestProjectionMappingTest:: testConsistentNonIncremental()
{
  preciceTrace ( "testConsistentNonIncremental()" );
  using namespace mesh;
  using utils::Vector2D;
  int dimensions = 2;

  // Create mesh to map from
  PtrMesh inMesh ( new Mesh("InMesh", dimensions, false) );
  PtrData inData = inMesh->createData ( "InData", 1 );
  int inDataID = inData->getID ();
  Vertex& v1 = inMesh->createVertex ( Vector2D(0.0, 0.0) );
  Vertex& v2 = inMesh->createVertex ( Vector2D(1.0, 1.0) );
  inMesh->createEdge ( v1, v2 );
  inMesh->computeState();
  inMesh->allocateDataValues();
  double valueVertex1 = 1.0;
  double valueVertex2 = 2.0;
  utils::DynVector& values = inData->values();
  values[0] = valueVertex1;
  values[1] = valueVertex2;

  // Create mesh to map to
  PtrMesh outMesh ( new Mesh("OutMesh", dimensions, false) );
  PtrData outData = outMesh->createData ( "OutData", 1 );
  int outDataID = outData->getID();

  // Setup mapping with mapping coordinates and geometry used
  NearestProjectionMapping mapping(Mapping::CONSISTENT, dimensions);
  mapping.setMeshes ( inMesh, outMesh );
  validateEquals ( mapping.hasComputedMapping(), false );

  Vertex& outv0 = outMesh->createVertex ( Vector2D(0.5, 0.5) );
  Vertex& outv1 = outMesh->createVertex ( Vector2D(-0.5, -0.5) );
  Vertex& outv2 = outMesh->createVertex ( Vector2D(1.5, 1.5) );
  outMesh->allocateDataValues();

  // Compute and perform mapping
  mapping.computeMapping();
  mapping.map ( inDataID, outDataID );

  // Validate results
  validateEquals ( mapping.hasComputedMapping(), true );
  validateNumericalEquals ( outData->values()[0], (valueVertex1 + valueVertex2) * 0.5 );
  validateNumericalEquals ( outData->values()[1], valueVertex1 );
  validateNumericalEquals ( outData->values()[2], valueVertex2 );

  // Redo mapping, results should be
  assign(outData->values()) = 0.0;
  mapping.map ( inDataID, outDataID );
  validateNumericalEquals ( outData->values()[0], (valueVertex1 + valueVertex2) * 0.5 );
  validateNumericalEquals ( outData->values()[1], valueVertex1 );
  validateNumericalEquals ( outData->values()[2], valueVertex2 );

  // Change vertex coordinates and redo mapping
  outv0.setCoords ( Vector2D(-0.5, -0.5) );
  outv1.setCoords ( Vector2D(1.5, 1.5) );
  outv2.setCoords ( Vector2D(0.5, 0.5) );
  assign(outData->values()) = 0.0;
  mapping.computeMapping();
  mapping.map ( inDataID, outDataID );
  validateNumericalEquals ( outData->values()[0], valueVertex1 );
  validateNumericalEquals ( outData->values()[1], valueVertex2 );
  validateNumericalEquals ( outData->values()[2], (valueVertex1 + valueVertex2) * 0.5 );

  // Reset output data to zero and redo the mapping
  assign(outData->values()) = 0.0;
  mapping.map ( inDataID, outDataID );
  validateNumericalEquals ( outData->values()[0], valueVertex1 );
  validateNumericalEquals ( outData->values()[1], valueVertex2 );
  validateNumericalEquals ( outData->values()[2], (valueVertex1 + valueVertex2) * 0.5 );
}

}}} // namespace precice, mapping, tests
