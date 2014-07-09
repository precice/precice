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
  if(doesFirstStep()){
    partnerNames.push_back(_secondParticipant);
  }
  else {
    partnerNames.push_back(_firstParticipant);
  }
  return partnerNames;
}


}} // namespace precice, cplscheme
