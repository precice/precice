// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "UncoupledScheme.hpp"
#include "com/Communication.hpp"

namespace precice {
namespace cplscheme {

tarch::logging::Log UncoupledScheme::
  _log ( "precice::cplscheme::UncoupledScheme" );

UncoupledScheme:: UncoupledScheme
(
  double             maxTime,
  int                maxTimesteps,
  int                validDigits,
  const std::string& participant )
:
  BaseCouplingScheme ( maxTime, maxTimesteps, UNDEFINED_TIMESTEP_LENGTH,
                   validDigits ),
  _participant ( participant )
{}

void UncoupledScheme:: initialize
(
  double startTime,
  int    startTimestep )
{
  preciceTrace2 ( "initialize()", startTime, startTimestep );
  setTime ( startTime );
  setTimesteps ( startTimestep );
  assertion1 ( tarch::la::greaterEquals(startTime, 0.0), startTime );
  assertion1 ( startTimestep >= 0, startTimestep );
  setIsInitialized(true);
}

void UncoupledScheme:: initializeData ()
{
  preciceTrace ( "initializeData()" );
}

void UncoupledScheme:: addComputedTime
(
  double timeToAdd )
{
  preciceTrace1("addComputedTime()", timeToAdd);
  preciceCheck ( isCouplingOngoing(), "addComputedTime()",
                 "Invalid call of addComputedTime() after simulation end!" );
# ifdef Asserts
  bool greaterThanZero = tarch::la::greater(timeToAdd, 0.0, _eps);
  assertion1(greaterThanZero, timeToAdd);
# endif // Asserts
  setComputedTimestepPart(getComputedTimestepPart() + timeToAdd);
  setTime(getTime() + timeToAdd);
}

void UncoupledScheme:: advance()
{
  preciceTrace("advance()");
  checkCompletenessRequiredActions();
  setIsCouplingTimestepComplete ( false );
  //double remainder = getTimestepRemainder ( computedTimestepLength );
  //setTime ( getTime() + computedTimestepLength );
  //if ( tarch::la::equals(getThisTimestepRemainder(), 0.0, eps) ){
  setIsCouplingTimestepComplete ( true );
  setTimesteps ( getTimesteps() + 1 );
    //if ( not tarch::la::equals(getMaxTime(),UNDEFINED_TIME) ){
    //  setMaxLengthNextTimestep ( getMaxTime() - getTime() );
    //}
  setComputedTimestepPart(0.0);
  //}
}

void UncoupledScheme:: finalize()
{
  preciceTrace ( "finalize()" );
}

std::vector<std::string> UncoupledScheme:: getCouplingPartners() const
{
  return std::vector<std::string>();
}

void UncoupledScheme:: sendState
(
  com::Communication::SharedPointer communication,
  int                   rankReceiver )
{
  communication->startSendPackage ( rankReceiver );
  BaseCouplingScheme::sendState ( communication, rankReceiver );
  communication->finishSendPackage();
}

void UncoupledScheme:: receiveState
(
  com::Communication::SharedPointer communication,
  int                   rankSender )
{
  communication->startReceivePackage ( rankSender );
  BaseCouplingScheme::receiveState ( communication, rankSender );
  communication->finishReceivePackage();
}

std::string UncoupledScheme:: printCouplingState () const
{
  std::ostringstream os;
  os << printBasicState() << " | " << printActionsState();
  return os.str();
}

}} // namespace precice, cplscheme
