// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
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
  preciceTrace ( "testMeasureData()" );

  // Create convergence measure for Vector data
  double convergenceLimit = 9.0;
  impl::AbsoluteConvergenceMeasure measure ( convergenceLimit );

  // Create data sets for old state of data and new state of data
  utils::DynVector oldValues0 ( 3 );
  utils::DynVector oldValues1 ( 3 );
  utils::DynVector oldValues2 ( 3 );
  utils::DynVector newValues  ( 3 );
  utils::DynVector designSpec ( 3, 0.0 );

  assignList(oldValues0) = -2.0, -1.0, 0.0;
  assignList(oldValues1) = 2.0, 3.0, 4.0;
  assignList(oldValues2) = 3.0, 4.0, 5.0;
  assignList(newValues) = 5.0, 6.0, 7.0;

  measure.measure ( oldValues0, newValues, designSpec );
  validate ( ! measure.isConvergence() );

  measure.measure ( oldValues1, newValues, designSpec );
  validate ( measure.isConvergence() );

  measure.measure ( oldValues2, newValues, designSpec );
  validate ( measure.isConvergence() );
}

}}} // namespace precice, cplscheme, tests
