#include "AbsoluteConvergenceMeasureTest.hpp"
#include "../impl/AbsoluteConvergenceMeasure.hpp"
#include "utils/Parallel.hpp"
#include "utils/Globals.hpp"

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::cplscheme::tests::AbsoluteConvergenceMeasureTest)

namespace precice {
namespace cplscheme {
namespace tests {

logging::Logger AbsoluteConvergenceMeasureTest::
   _log ( "precice::cplscheme::tests::AbsoluteConvergenceMeasureTest" );

AbsoluteConvergenceMeasureTest:: AbsoluteConvergenceMeasureTest ()
:
  TestCase ( "cplscheme::tests::AbsoluteConvergenceMeasureTest" )
{}

void AbsoluteConvergenceMeasureTest:: run ()
{
  PRECICE_MASTER_ONLY {
    testMethod ( testMeasureData );
  }
}

void AbsoluteConvergenceMeasureTest:: testMeasureData ()
{
  TRACE();
  using Eigen::Vector3d;
  // Create convergence measure for Vector data
  double convergenceLimit = 9.0;
  impl::AbsoluteConvergenceMeasure measure ( convergenceLimit );

  // Create data sets for old state of data and new state of data
  Vector3d oldValues0(-2, -1, 0);
  Vector3d oldValues1(2, 3, 4);
  Vector3d oldValues2(3, 4, 5);
  Vector3d newValues(5, 6, 7);
  Vector3d designSpec = Vector3d::Zero();

  measure.measure ( oldValues0, newValues, designSpec );
  validate ( ! measure.isConvergence() );

  measure.measure ( oldValues1, newValues, designSpec );
  validate ( measure.isConvergence() );

  measure.measure ( oldValues2, newValues, designSpec );
  validate ( measure.isConvergence() );
}

}}} // namespace precice, cplscheme, tests
