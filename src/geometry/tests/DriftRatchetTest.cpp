#include "DriftRatchetTest.hpp"
#include "geometry/DriftRatchet.hpp"
#include "mesh/Mesh.hpp"
#include "io/ExportVTK.hpp"
#include "mesh/SharedPointer.hpp"
#include "utils/Parallel.hpp"
#include "utils/Globals.hpp"
#include "boost/smart_ptr.hpp"

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::geometry::tests::DriftRatchetTest)

namespace precice {
namespace geometry {
namespace tests {

logging::Logger DriftRatchetTest:: _log ( "precice::geometry::tests::DriftRatchetTest" );

DriftRatchetTest:: DriftRatchetTest (void)
:
  TestCase ("geometry::DriftRatchetTest")
{}

void DriftRatchetTest:: run ()
{
  PRECICE_MASTER_ONLY {
    preciceTrace ( "run" );
    for ( int dim=2; dim <= 3; dim++ ){
      bool flipNormals = true;
      mesh::Mesh mesh ( "test-driftratchet", dim, flipNormals );
      const double pores = 2.0;
      double maxRadius = 2.4;
      double length = geometry::DriftRatchet::getCharacteristicLength(maxRadius) * pores;
      double discretizationWidth = 0.5;
      double minRadius = geometry::DriftRatchet::getDefaultMinRadius(maxRadius);
      geometry::DriftRatchet driftRatchet (
          utils::DynVector(dim,0.0), discretizationWidth, maxRadius, minRadius,
          geometry::DriftRatchet::getDefaultShapeParameter(),
          length, pores, 0, 1, 2 );
      driftRatchet.create ( mesh );
      preciceDebug ( "Created Container with " << mesh.triangles().size()
                     << " triangles and " << mesh.vertices().size() << " vertices" );
      io::ExportVTK exportVTK(true);
      std::ostringstream filename;
      filename << "geometry-DriftRatchetTest-" << dim;
      exportVTK.doExport ( filename.str(), mesh );
    }
  }
}

}}} // namespace precice, geometry, tests
