#include "SphereTest.hpp"
#include "geometry/Sphere.hpp"
#include "io/ExportVTK.hpp"
#include "mesh/Mesh.hpp"
#include "utils/Parallel.hpp"
#include <sstream>

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::geometry::tests::SphereTest)

namespace precice {
namespace geometry {
namespace tests {

logging::Logger SphereTest:: _log ( "precice::geometry::tests::SphereTest" );

SphereTest:: SphereTest ()
:
  TestCase ("SphereTest")
{}

void SphereTest:: run ()
{
  PRECICE_MASTER_ONLY {
    preciceTrace ( "run" );
    for ( int dim=2; dim <= 3; dim++ ){
      bool flipNormals = false;
      mesh::Mesh mesh ( "test-sphere", dim, flipNormals );
      Eigen::VectorXd offset = Eigen::VectorXd::Zero(dim);
      double discretizationWidth = 0.1;
      double radius = 1.0;
      Sphere sphere ( offset, discretizationWidth, radius );
      sphere.create ( mesh );
      io::ExportVTK exportVTK(true);
      std::ostringstream filename;
      filename << "geometry-SphereTest-" << dim;
      std::string location = "";
      exportVTK.doExport ( filename.str(), location, mesh );
    }
  }
}

}}} // namespace precice, geometry, tests
