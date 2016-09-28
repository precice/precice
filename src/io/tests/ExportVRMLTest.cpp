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

logging::Logger ExportVRMLTest:: _log ("precice::io::ExportVRMLTest");

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
    std::string location = "";
    ex.doExport( stream.str(), location, mesh );
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

    std::string name ( "side-0" );
    mesh.setSubID ( name );

    name = "side-1";
    mesh.setSubID ( name );

    name = "side-2";
    mesh.setSubID ( name );

    name = "side-3";
    mesh.setSubID ( name );
    mesh.createData ( "Data", dim );
    cuboid.create ( mesh );

    DEBUG ( "Mesh vertices = " << mesh.vertices().size() );
    std::ostringstream stream;
    stream << "io-ExportVRMLTest-testExportCuboid-" << dim << "d.wrl";
    ExportVRML ex(false);
    std::string location = "";
    ex.doExport ( stream.str(), location, mesh );
  }
}

}}} // namespace precice, io, tests
