// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "UncoupledCouplingScheme.hpp"
#include "com/Communication.hpp"
#include <limits>

namespace precice {
namespace cplscheme {

tarch::logging::Log UncoupledCouplingScheme::
  _log ( "precice::cplscheme::UncoupledCouplingScheme" );

UncoupledCouplingScheme:: UncoupledCouplingScheme
(
  double             maxTime,
  int                maxTimesteps,
  int                validDigits,
  const std::string& participant )
:
  CouplingScheme ( maxTime, maxTimesteps, UNDEFINED_TIMESTEP_LENGTH,
                   validDigits ),
  _participant ( participant )
{}

void UncoupledCouplingScheme:: initialize
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

void UncoupledCouplingScheme:: initializeData ()
{
  preciceTrace ( "initializeData()" );
}

void UncoupledCouplingScheme:: addComputedTime
(
  double timeToAdd )
{
  preciceTrace1("addComputedTime()", timeToAdd);
  preciceCheck ( isCouplingOngoing(), "addComputedTime()",
                 "Invalid call of addComputedTime() after simulation end!" );
# ifdef Asserts
  double eps = std::pow ( 10.0, -1 * getValidDigits() );
  bool greaterThanZero = tarch::la::greater(timeToAdd, 0.0, eps);
  assertion1(greaterThanZero, timeToAdd);
# endif // Asserts
  setComputedTimestepPart(getComputedTimestepPart() + timeToAdd);
  setTime(getTime() + timeToAdd);
}

void UncoupledCouplingScheme:: advance()
{
  preciceTrace("advance()");
  checkCompletenessRequiredActions();
  setIsCouplingTimestepComplete ( false );
  //double remainder = getTimestepRemainder ( computedTimestepLength );
  //setTime ( getTime() + computedTimestepLength );
  //double eps = std::pow ( 10.0, -1 * getValidDigits() );
  //if ( tarch::la::equals(getThisTimestepRemainder(), 0.0, eps) ){
  setIsCouplingTimestepComplete ( true );
  setTimesteps ( getTimesteps() + 1 );
    //if ( not tarch::la::equals(getMaxTime(),UNDEFINED_TIME) ){
    //  setMaxLengthNextTimestep ( getMaxTime() - getTime() );
    //}
  setComputedTimestepPart(0.0);
  //}
}

void UncoupledCouplingScheme:: finalize()
{
  preciceTrace ( "finalize()" );
}

std::vector<std::string> UncoupledCouplingScheme:: getCouplingPartners() const
{
  return std::vector<std::string>();
}

void UncoupledCouplingScheme:: sendState
(
  com::PtrCommunication communication,
  int                   rankReceiver )
{
  communication->startSendPackage ( rankReceiver );
  CouplingScheme::sendState ( communication, rankReceiver );
  communication->finishSendPackage();
}

void UncoupledCouplingScheme:: receiveState
(
  com::PtrCommunication communication,
  int                   rankSender )
{
  communication->startReceivePackage ( rankSender );
  CouplingScheme::receiveState ( communication, rankSender );
  communication->finishReceivePackage();
}

std::string UncoupledCouplingScheme:: printCouplingState () const
{
  std::ostringstream os;
  os << printBasicState() << " | " << printActionsState();
  return os.str();
}

}} // namespace precice, cplscheme
