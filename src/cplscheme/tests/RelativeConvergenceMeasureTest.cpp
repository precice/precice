#include "RelativeConvergenceMeasureTest.hpp"
#include "../impl/RelativeConvergenceMeasure.hpp"
#include "utils/Parallel.hpp"
#include "utils/Globals.hpp"

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::cplscheme::tests::RelativeConvergenceMeasureTest)

namespace precice {
namespace cplscheme {
namespace tests {

logging::Logger RelativeConvergenceMeasureTest::
   _log ( "precice::cplscheme::tests::RelativeConvergenceMeasureTest" );

RelativeConvergenceMeasureTest:: RelativeConvergenceMeasureTest ()
:
  TestCase ( "precice::cplscheme::tests::RelativeConvergenceMeasureTest" )
{}

void RelativeConvergenceMeasureTest:: run ()
{
  PRECICE_MASTER_ONLY {
    testMethod ( testMeasureData );
  }
}

void RelativeConvergenceMeasureTest:: testMeasureData ()
{
  TRACE();
  using Eigen::Vector3d;
  double convergenceLimit = 0.1; // 10%
  impl::RelativeConvergenceMeasure measure ( convergenceLimit );

  // Create data sets for old state of data and new state of data
  Vector3d oldValues0 (1, 1, 1);
  Vector3d oldValues1 (2, 2, 2);
  Vector3d oldValues2 (2.9, 2.9, 2.9);
  Vector3d newValues (3, 3, 3);
  Vector3d designSpec ( 0, 0, 0 );
  
  measure.measure ( oldValues0, newValues, designSpec );
  validate ( ! measure.isConvergence() );

  measure.measure ( oldValues1, newValues, designSpec );
  validate ( ! measure.isConvergence() );

  measure.measure ( oldValues2, newValues, designSpec );
  validate ( measure.isConvergence() );
}

}}} // namespace precice, cplscheme, tests
