// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "ExplicitCouplingScheme.hpp"
#include "com/Communication.hpp"
#include "utils/Globals.hpp"
#include "utils/Parallel.hpp"

namespace precice {
namespace cplscheme {

tarch::logging::Log ExplicitCouplingScheme::
   _log("precice::cplscheme::ExplicitCouplingScheme");

ExplicitCouplingScheme:: ExplicitCouplingScheme
(
  double                maxTime,
  int                   maxTimesteps,
  double                timestepLength,
  int                   validDigits,
  const std::string&    firstParticipant,
  const std::string&    secondParticipant,
  const std::string&    localParticipantName,
  com::PtrCommunication communication,
  constants::TimesteppingMethod dtMethod)
:
  CouplingScheme(maxTime, maxTimesteps, timestepLength, validDigits),
  _firstParticipant(firstParticipant),
  _secondParticipant(secondParticipant),
  _doesFirstStep(false),
  _communication(communication),
  _participantSetsDt(false),
  _participantReceivesDt(false)
{
  preciceCheck(_firstParticipant != _secondParticipant,
               "ExplicitCouplingScheme()", "First participant and "
               << "second participant must be different!");
  if (dtMethod == constants::FIXED_DT){
    preciceCheck(hasTimestepLength(), "ExplicitCouplingScheme()",
        "Timestep length value has to be given "
        << "when the fixed timestep length method is chosen for an explicit "
        << "coupling scheme!");
  }
  if (localParticipantName == _firstParticipant){
    _doesFirstStep = true;
    if (dtMethod == constants::FIRST_PARTICIPANT_SETS_DT){
      _participantSetsDt = true;
      setTimestepLength(UNDEFINED_TIMESTEP_LENGTH);
    }
  }
  else if (localParticipantName == _secondParticipant){
    if (dtMethod == constants::FIRST_PARTICIPANT_SETS_DT){
      _participantReceivesDt = true;
    }
  }
  else {
    preciceError("initialize()", "Name of local participant \""
                 << localParticipantName << "\" does not match any "
                 << "participant specified for the coupling scheme!");
  }
  assertion(communication.use_count() > 0);
}

ExplicitCouplingScheme:: ~ExplicitCouplingScheme()
{
}

void ExplicitCouplingScheme:: initialize
(
  double startTime,
  int    startTimestep)
{
  preciceTrace2("initialize()", startTime, startTimestep);
  assertion1(tarch::la::greaterEquals(startTime, 0.0), startTime);
  assertion1(startTimestep >= 0, startTimestep);
  assertion(_communication->isConnected());
  setTime(startTime);
  setTimesteps(startTimestep);
  if((not _doesFirstStep) && isCouplingOngoing()){
    preciceDebug("Receiving data...");
    _communication->startReceivePackage(0);
    if (_participantReceivesDt){
      double dt = UNDEFINED_TIMESTEP_LENGTH;
      _communication->receive(dt, 0);
      assertion(not tarch::la::equals(dt, UNDEFINED_TIMESTEP_LENGTH));
      setTimestepLength(dt);
    }
    receiveData(_communication);
    _communication->finishReceivePackage();
    setHasDataBeenExchanged(true);
  }
  setIsInitialized(true);
}

void ExplicitCouplingScheme:: addComputedTime
(
  double timeToAdd)
{
  preciceTrace1("addComputedTime()", timeToAdd);
  preciceCheck(isCouplingOngoing(), "addComputedTime()",
                 "Invalid call of addComputedTime() after simulation end!");

  // Check validness
  double eps = std::pow(10.0, -1 * getValidDigits());
  bool greaterThanZero = tarch::la::greater(timeToAdd, 0.0, eps);
  preciceCheck(greaterThanZero, "addComputedTime()", "The computed timestep length "
               << "exceeds the maximum timestep limit for this time step!");

  setComputedTimestepPart(getComputedTimestepPart() + timeToAdd);
  setTime(getTime() + timeToAdd);
}

void ExplicitCouplingScheme:: advance()
{
  preciceTrace("advance()");
  checkCompletenessRequiredActions();
  setHasDataBeenExchanged(false);
  setIsCouplingTimestepComplete(false);
  //double remainder = getTimestepRemainder();
  //setTime(getTime() + computedTimestepLength);
  double eps = std::pow(10.0, -1 * getValidDigits());
  if (tarch::la::equals(getThisTimestepRemainder(), 0.0, eps)){
    setIsCouplingTimestepComplete(true);
    setTimesteps(getTimesteps() + 1);
    //if (isCouplingOngoing() || _doesFirstStep){
    preciceDebug("Sending data...");
    _communication->startSendPackage(0);
    if (_participantSetsDt){
      _communication->send(getComputedTimestepPart(), 0);
    }
    sendData(_communication);
    _communication->finishSendPackage();
    //}
    if (isCouplingOngoing() || _doesFirstStep){
      preciceDebug("Receiving data...");
      _communication->startReceivePackage(0);
      if (_participantReceivesDt){
        double dt = UNDEFINED_TIMESTEP_LENGTH;
        _communication->receive(dt, 0);
        assertion(not tarch::la::equals(dt, UNDEFINED_TIMESTEP_LENGTH));
        setTimestepLength(dt);
      }
      receiveData(_communication);
      _communication->finishReceivePackage();
    }
    setHasDataBeenExchanged(true);
    setComputedTimestepPart(0.0);
    //setMaxLengthNextTimestep(getTimestepLength());
  }
}

void ExplicitCouplingScheme:: finalize()
{
   preciceTrace("finalize()");
   checkCompletenessRequiredActions();
}

void ExplicitCouplingScheme:: sendState
(
  com::PtrCommunication communication,
  int                   rankReceiver)
{
  communication->startSendPackage(0);
  CouplingScheme::sendState(communication, rankReceiver);
  communication->finishSendPackage();
}

void ExplicitCouplingScheme:: receiveState
(
  com::PtrCommunication communication,
  int                   rankSender)
{
  communication->startSendPackage(0);
  CouplingScheme::receiveState(communication, rankSender);
  communication->finishSendPackage();
}

std::string ExplicitCouplingScheme:: printCouplingState() const
{
   std::ostringstream os;
   os << printBasicState() << " | " << printActionsState();
   return os.str();
}

std::vector<std::string> ExplicitCouplingScheme:: getCouplingPartners
(
  const std::string& accessorName) const
{
  std::vector<std::string> partnerNames;
  if(accessorName == _firstParticipant){
    partnerNames.push_back(_secondParticipant);
  }
  else if(accessorName == _secondParticipant){
    partnerNames.push_back(_firstParticipant);
  }
  else {
    preciceError("getCouplingPartners()", "No coupling partner could be found.");
  }
  return partnerNames;
}


}} // namespace precice, cplscheme
