// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "UncoupledCouplingSchemeTest.hpp"
#include "cplscheme/UncoupledCouplingScheme.hpp"
#include "utils/Parallel.hpp"

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::cplscheme::tests::UncoupledCouplingSchemeTest)

namespace precice {
namespace cplscheme {
namespace tests {

tarch::logging::Log UncoupledCouplingSchemeTest::
   _log ( "precice::cplscheme::tests::UncoupledCouplingSchemeTest" );

UncoupledCouplingSchemeTest:: UncoupledCouplingSchemeTest ()
:
   TestCase ( "precice::cplscheme::tests::UncoupledCouplingSchemeTest" )
{}

void UncoupledCouplingSchemeTest:: run ()
{
  PRECICE_MASTER_ONLY {
    testMethod ( testBasics );
  }
}

void UncoupledCouplingSchemeTest:: testBasics ()
{
  preciceTrace ( "testBasics()" );
  typedef UncoupledCouplingScheme CplScheme;
  double maxTime = 10.0;
  int maxTimesteps = CplScheme::UNDEFINED_TIMESTEPS;
  UncoupledCouplingScheme cplscheme ( maxTime, maxTimesteps, 14, "TestParticipant" );
  cplscheme.initialize ( 0.0, 0 );
  validate ( not cplscheme.isCouplingTimestepComplete() );

  cplscheme.addComputedTime(1.0);
  cplscheme.advance();
  using tarch::la::equals;
  double maxDt = cplscheme.getNextTimestepMaxLength();
  validateWithParams1 ( equals(maxDt, 9.0), maxDt );
  validate ( cplscheme.isCouplingTimestepComplete() );

  cplscheme.addComputedTime(2.0);
  maxDt = cplscheme.getNextTimestepMaxLength();
  validateWithParams1 ( equals(maxDt, 7.0), maxDt );
  validate ( cplscheme.isCouplingTimestepComplete() );

  cplscheme.addComputedTime(3.0);
  cplscheme.advance();
  maxDt = cplscheme.getNextTimestepMaxLength();
  validateWithParams1 ( equals(maxDt, 4.0), maxDt );

  cplscheme.addComputedTime(4.0);
  cplscheme.advance();
  maxDt = cplscheme.getNextTimestepMaxLength();
  validateWithParams1 ( equals(maxDt, 0.0), maxDt );
  validate ( not cplscheme.isCouplingOngoing() );
  validate ( cplscheme.isCouplingTimestepComplete() );

  cplscheme.finalize();
}

}}} // namespace precice, cplscheme, tests
