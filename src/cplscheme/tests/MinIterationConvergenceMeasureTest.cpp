// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
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

tarch::logging::Log MinIterationConvergenceMeasureTest::
   _log ( "precice::cplscheme::tests::MinIterationConvergenceMeasureTest" );


MinIterationConvergenceMeasureTest:: MinIterationConvergenceMeasureTest ()
:
  TestCase ( "precice::cplscheme::tests::MinIterationConvergenceMeasureTest" )
{}

void MinIterationConvergenceMeasureTest:: run ()
{
  PRECICE_MASTER_ONLY {
    preciceTrace ( "run()" );
    impl::MinIterationConvergenceMeasure measure ( 5 );
    utils::DynVector emptyValues; // No values needed for min-iter

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
