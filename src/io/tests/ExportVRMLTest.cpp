#include "ExportVRMLTest.hpp"
#include "io/ExportVRML.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "utils/Parallel.hpp"
#include <string>

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::io::tests::ExportVRMLTest)

namespace precice {
namespace io {
namespace tests {

using namespace std;

logging::Logger ExportVRMLTest:: _log ("io::ExportVRMLTest");

ExportVRMLTest:: ExportVRMLTest()
:
   TestCase ("io::ExportVRMLTest")
{}

void ExportVRMLTest:: run()
{
  PRECICE_MASTER_ONLY {
    testMethod ( testExportSimpleMesh );
  }
}


void ExportVRMLTest:: testExportSimpleMesh()
{
  TRACE();
  int dim=2;
  mesh::Mesh::resetGeometryIDsGlobally ();
  bool flipNormals = false;
  mesh::Mesh mesh ( "test-cuboid", dim, flipNormals );

  // some dummy mesh
  mesh::Vertex& v1 = mesh.createVertex(Eigen::Vector2d(1.0, 1.0));
  mesh::Vertex& v2 = mesh.createVertex(Eigen::Vector2d(2.0, 1.0));
  mesh::Vertex& v3 = mesh.createVertex(Eigen::Vector2d(1.0, 2.0));
  mesh::Vertex& v4 = mesh.createVertex(Eigen::Vector2d(2.0, 2.0));
  mesh.createEdge(v1,v2);
  mesh.createEdge(v1,v3);
  mesh.createEdge(v2,v4);
  mesh.createEdge(v3,v4);

  mesh.createData ( "Data", dim );
  mesh.allocateDataValues();

  DEBUG ( "Mesh vertices = " << mesh.vertices().size() );
  std::ostringstream stream;
  stream << "io-ExportVRMLTest-testExportCuboid-" << dim << "d.wrl";
  ExportVRML ex(false);
  std::string location = "";
  ex.doExport ( stream.str(), location, mesh );
}

}}} // namespace precice, io, tests
