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
  BaseCouplingScheme(maxTime, maxTimesteps, timestepLength, validDigits),
  _firstParticipant(firstParticipant),
  _secondParticipant(secondParticipant),
  _doesFirstStep(false),
  _communication(communication),
  _participantSetsDt(false),
  _participantReceivesDt(false),
  _hasToReceiveInitData(false),
  _hasToSendInitData(false)
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
  assertion(not isInitialized());
  assertion1(tarch::la::greaterEquals(startTime, 0.0), startTime);
  assertion1(startTimestep >= 0, startTimestep);
  assertion(getCommunication()->isConnected());
  preciceCheck(not getSendData().empty(), "initialize()",
               "No send data configured! Use explicit scheme for one-way coupling.");
  setTime(startTime);
  setTimesteps(startTimestep);


  foreach (DataMap::value_type & pair, getSendData()){
    if (pair.second->initialize){
      preciceCheck(not doesFirstStep(), "initialize()",
                   "Only second participant can initialize data!");
      preciceDebug("Initialized data to be written");
      setHasToSendInitData(true);
      break;
    }
  }

  foreach (DataMap::value_type & pair, getReceiveData()){
    if (pair.second->initialize){
      preciceCheck(doesFirstStep(), "initialize()",
                   "Only first participant can receive initial data!");
      preciceDebug("Initialized data to be received");
      setHasToReceiveInitData(true);
    }
  }

   // If the second participant initializes data, the first receive for the
   // second participant is done in initializeData() instead of initialize().
  if ((not doesFirstStep()) && (not hasToSendInitData()) && isCouplingOngoing()){
    preciceDebug("Receiving data");
    getCommunication()->startReceivePackage(0);
    if (participantReceivesDt()){
      double dt = UNDEFINED_TIMESTEP_LENGTH;
      getCommunication()->receive(dt, 0);
      preciceDebug("received timestep length of " << dt);
      assertion(not tarch::la::equals(dt, UNDEFINED_TIMESTEP_LENGTH));
      setTimestepLength(dt);
    }
    receiveData(getCommunication());
    getCommunication()->finishReceivePackage();
    setHasDataBeenExchanged(true);
  }

  if(hasToSendInitData()){
    requireAction(constants::actionWriteInitialData());
  }
  setIsInitialized(true);
}

void ExplicitCouplingScheme:: initializeData()
{
  preciceTrace("initializeData()");
  preciceCheck(isInitialized(), "initializeData()",
     "initializeData() can be called after initialize() only!");

  if((not hasToSendInitData()) && (not hasToReceiveInitData())){
    preciceInfo("initializeData()", "initializeData is skipped since no data has to be initialized");
    return;
  }

  preciceDebug("Initializing Data ...");

  preciceCheck(not (hasToSendInitData() && isActionRequired(constants::actionWriteInitialData())),
     "initializeData()", "InitialData has to be written to preCICE before calling initializeData()");

  setHasDataBeenExchanged(false);

  if (hasToReceiveInitData() && isCouplingOngoing()){
    assertion(doesFirstStep());
    preciceDebug("Receiving data");
    getCommunication()->startReceivePackage(0);
    if (participantReceivesDt()){
      double dt = UNDEFINED_TIMESTEP_LENGTH;
      getCommunication()->receive(dt, 0);
      preciceDebug("received timestep length of " << dt);
      assertion(not tarch::la::equals(dt, UNDEFINED_TIMESTEP_LENGTH));
      setTimestepLength(dt);
      //setMaxLengthNextTimestep(dt);
    }
    receiveData(getCommunication());
    getCommunication()->finishReceivePackage();
    setHasDataBeenExchanged(true);
  }



  if (hasToSendInitData() && isCouplingOngoing()){
    assertion(not doesFirstStep());

    // The second participant sends the initialized data to the first particpant
    // here, which receives the data on call of initialize().
    sendData(getCommunication());
    getCommunication()->startReceivePackage(0);
    // This receive replaces the receive in initialize().
    receiveData(getCommunication());
    getCommunication()->finishReceivePackage();
    setHasDataBeenExchanged(true);

  }

  //in order to check in advance if initializeData has been called (if necessary)
  setHasToSendInitData(false);
  setHasToReceiveInitData(false);
}

void ExplicitCouplingScheme:: advance()
{
  preciceTrace("advance()");
  checkCompletenessRequiredActions();
  setHasDataBeenExchanged(false);
  setIsCouplingTimestepComplete(false);
  double eps = std::pow(10.0, -1 * getValidDigits());
  if (tarch::la::equals(getThisTimestepRemainder(), 0.0, eps)){
    setIsCouplingTimestepComplete(true);
    setTimesteps(getTimesteps() + 1);
    preciceDebug("Sending data...");
    getCommunication()->startSendPackage(0);
    if (participantSetsDt()){
      getCommunication()->send(getComputedTimestepPart(), 0);
    }
    sendData(getCommunication());
    getCommunication()->finishSendPackage();

    if (isCouplingOngoing() || doesFirstStep()){
      preciceDebug("Receiving data...");
      getCommunication()->startReceivePackage(0);
      if (participantReceivesDt()){
        double dt = UNDEFINED_TIMESTEP_LENGTH;
        getCommunication()->receive(dt, 0);
        assertion(not tarch::la::equals(dt, UNDEFINED_TIMESTEP_LENGTH));
        setTimestepLength(dt);
      }
      receiveData(getCommunication());
      getCommunication()->finishReceivePackage();
    }
    setHasDataBeenExchanged(true);
    setComputedTimestepPart(0.0);
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
  BaseCouplingScheme::sendState(communication, rankReceiver);
  communication->finishSendPackage();
}

void ExplicitCouplingScheme:: receiveState
(
  com::PtrCommunication communication,
  int                   rankSender)
{
  communication->startSendPackage(0);
  BaseCouplingScheme::receiveState(communication, rankSender);
  communication->finishSendPackage();
}

std::string ExplicitCouplingScheme:: printCouplingState() const
{
   std::ostringstream os;
   os << printBasicState() << " | " << printActionsState();
   return os.str();
}

std::vector<std::string> ExplicitCouplingScheme:: getCouplingPartners() const
{
  std::vector<std::string> partnerNames;
  // Add non-local participant
  if(_doesFirstStep){
    partnerNames.push_back(_secondParticipant);
  }
  else {
    partnerNames.push_back(_firstParticipant);
  }
  return partnerNames;
}


}} // namespace precice, cplscheme
