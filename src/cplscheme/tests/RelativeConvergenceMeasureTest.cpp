// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "RelativeConvergenceMeasureTest.hpp"
#include "../impl/RelativeConvergenceMeasure.hpp"
#include "utils/Parallel.hpp"
#include "utils/Globals.hpp"

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::cplscheme::tests::RelativeConvergenceMeasureTest)

namespace precice {
namespace cplscheme {
namespace tests {

tarch::logging::Log RelativeConvergenceMeasureTest::
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
  preciceTrace ( "testMeasureVectorData()" );
  double convergenceLimit = 0.1; // 10%
  impl::RelativeConvergenceMeasure measure ( convergenceLimit );

  // Create data sets for old state of data and new state of data
  utils::DynVector oldValues0 ( 3 );
  utils::DynVector oldValues1 ( 3 );
  utils::DynVector oldValues2 ( 3 );
  utils::DynVector newValues ( 3 );
  utils::DynVector designSpec ( 3, 0.0 );
  assignList(oldValues0) = 1.0, 1.0, 1.0;
  assignList(oldValues1) = 2.0, 2.0, 2.0;
  assignList(oldValues2) = 2.9, 2.9, 2.9;
  assignList(newValues) = 3.0, 3.0, 3.0;

  measure.measure ( oldValues0, newValues, designSpec );
  validate ( ! measure.isConvergence() );

  measure.measure ( oldValues1, newValues, designSpec );
  validate ( ! measure.isConvergence() );

  measure.measure ( oldValues2, newValues, designSpec );
  validate ( measure.isConvergence() );
}

}}} // namespace precice, cplscheme, tests
