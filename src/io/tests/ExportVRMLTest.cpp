// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "ExportVRMLTest.hpp"
#include "io/ExportVRML.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Data.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Mesh.hpp"
#include "geometry/DriftRatchet.hpp"
#include "geometry/Cuboid.hpp"
#include "utils/Parallel.hpp"
#include "utils/Dimensions.hpp"
#include "utils/Globals.hpp"
#include <iostream>
#include <map>
#include <string>

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::io::tests::ExportVRMLTest)

namespace precice {
namespace io {
namespace tests {

using namespace std;

tarch::logging::Log ExportVRMLTest:: _log ("precice::io::ExportVRMLTest");

ExportVRMLTest:: ExportVRMLTest()
:
   TestCase ("io::ExportVRMLTest")
{}

void ExportVRMLTest:: run()
{
  PRECICE_MASTER_ONLY {
    testMethod ( testExportDriftRatchet );
    testMethod ( testExportCuboid );
  }
}

void ExportVRMLTest:: testExportDriftRatchet()
{
  preciceTrace ( "testExportDriftRatchet" );
  for ( int dim=2; dim <= 3; dim++ ){
    ExportVRML ex(false);
    bool flipNormals = false;
    mesh::Mesh mesh ( "Ratchet", dim, flipNormals );
    geometry::DriftRatchet ratchet (
      utils::DynVector(dim, 0.0), // offset
      0.1,         // discretization width
      geometry::DriftRatchet::getDefaultMaxRadius(),
      geometry::DriftRatchet::getDefaultMinRadius(geometry::DriftRatchet::getDefaultMaxRadius()),
      geometry::DriftRatchet::getDefaultShapeParameter(),
      geometry::DriftRatchet::getCharacteristicLength(geometry::DriftRatchet::getDefaultMaxRadius()),
      1.5, // pores
      0,   // wall index
      1,   // inflow index
      2 ); // outflow index
    mesh::PtrData data = mesh.createData ( "Data", dim );
    ratchet.create ( mesh );
    std::ostringstream stream;
    stream << "io-ExportVRMLTest-testExportDriftRatchet-" << dim << "d.wrl";
    ex.doExport( stream.str(), mesh );
  }
}

void ExportVRMLTest:: testExportCuboid()
{
  preciceTrace ( "testExportCuboid" );
  for ( int dim=2; dim <= 3; dim++ ){
    mesh::Mesh::resetGeometryIDsGlobally ();
    bool flipNormals = false;
    mesh::Mesh mesh ( "test-cuboid", dim, flipNormals );
    utils::DynVector offset(dim, 0.0);
    double h = 1.0;
    utils::DynVector length(dim, 5.0);
    geometry::Cuboid cuboid ( offset, h, length );

    int sideIndex = 0;
    std::string name ( "side-0" );
    mesh.setSubID ( name );

    sideIndex = 1;
    name = "side-1";
    mesh.setSubID ( name );

    sideIndex = 2;
    name = "side-2";
    mesh.setSubID ( name );

    sideIndex = 3;
    name = "side-3";
    mesh.setSubID ( name );
    mesh.createData ( "Data", dim );
    cuboid.create ( mesh );

    preciceDebug ( "Mesh vertices = " << mesh.vertices().size() );
    std::ostringstream stream;
    stream << "io-ExportVRMLTest-testExportCuboid-" << dim << "d.wrl";
    ExportVRML ex(false);
    ex.doExport ( stream.str(), mesh );
  }
}

}}} // namespace precice, io, tests
