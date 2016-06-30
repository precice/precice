#include "UncoupledSchemeTest.hpp"
#include "cplscheme/UncoupledScheme.hpp"
#include "utils/Parallel.hpp"

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::cplscheme::tests::UncoupledSchemeTest)

namespace precice {
namespace cplscheme {
namespace tests {

logging::Logger UncoupledSchemeTest::
   _log ( "precice::cplscheme::tests::UncoupledSchemeTest" );

UncoupledSchemeTest:: UncoupledSchemeTest ()
:
   TestCase ( "precice::cplscheme::tests::UncoupledSchemeTest" )
{}

void UncoupledSchemeTest:: run ()
{
  PRECICE_MASTER_ONLY {
    testMethod ( testBasics );
  }
}

void UncoupledSchemeTest:: testBasics ()
{
  preciceTrace ( "testBasics()" );
  typedef UncoupledScheme CplScheme;
  double maxTime = 10.0;
  int maxTimesteps = CplScheme::UNDEFINED_TIMESTEPS;
  UncoupledScheme cplscheme ( maxTime, maxTimesteps, 14, "TestParticipant" );
  cplscheme.initialize ( 0.0, 1 );
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
