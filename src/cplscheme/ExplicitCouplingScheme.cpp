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
  BaseCouplingScheme(maxTime, maxTimesteps, timestepLength, validDigits,
		     firstParticipant, secondParticipant, localParticipantName,
		     communication, 1, dtMethod)
{}


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
