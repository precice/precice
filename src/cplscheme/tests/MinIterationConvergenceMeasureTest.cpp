#include "MinIterationConvergenceMeasureTest.hpp"
#include "../impl/MinIterationConvergenceMeasure.hpp"
#include "utils/Globals.hpp"
#include "utils/Parallel.hpp"
#include "mesh/Mesh.hpp"

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::cplscheme::tests::MinIterationConvergenceMeasureTest)

namespace precice {
namespace cplscheme {
namespace tests {

logging::Logger MinIterationConvergenceMeasureTest::
   _log ( "precice::cplscheme::tests::MinIterationConvergenceMeasureTest" );


MinIterationConvergenceMeasureTest:: MinIterationConvergenceMeasureTest ()
:
  TestCase ( "precice::cplscheme::tests::MinIterationConvergenceMeasureTest" )
{}

void MinIterationConvergenceMeasureTest:: run ()
{
  PRECICE_MASTER_ONLY {
    TRACE();
    impl::MinIterationConvergenceMeasure measure ( 5 );
    Eigen::VectorXd emptyValues; // No values needed for min-iter

    for ( int iSeries=0; iSeries < 3; iSeries++ ) {
      measure.newMeasurementSeries ();
      for ( int iMeasurement=1; iMeasurement < 10; iMeasurement++ ) {
        measure.measure ( emptyValues, emptyValues, emptyValues );
        if ( iMeasurement < 5 ) {
          validate ( ! measure.isConvergence() );
        }
        else {
          validate ( measure.isConvergence() );
        }
      }
    }
  }
}

}}} // namespace precice, cplscheme, tests
