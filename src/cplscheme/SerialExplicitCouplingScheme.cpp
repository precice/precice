// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "SerialExplicitCouplingScheme.hpp"
#include "Constants.hpp"
#include "com/Communication.hpp"
#include "com/SharedPointer.hpp"

namespace precice {
namespace cplscheme {

tarch::logging::Log SerialExplicitCouplingScheme::
_log("precice::cplscheme::SerialExplicitCouplingScheme" );

SerialExplicitCouplingScheme:: SerialExplicitCouplingScheme
(
  double                maxTime,
  int                   maxTimesteps,
  double                timestepLength,
  int                   validDigits,
  const std::string&    firstParticipant,
  const std::string&    secondParticipant,
  const std::string&    localParticipant,
  com::PtrCommunication communication,
  constants::TimesteppingMethod dtMethod )
  :
  ExplicitCouplingScheme(maxTime,maxTimesteps,timestepLength,validDigits,firstParticipant,
			 secondParticipant,localParticipant,communication,dtMethod)
{}

void SerialExplicitCouplingScheme:: advance()
{
  preciceTrace2("advance()", getTimesteps(), getTime());
  checkCompletenessRequiredActions();
  preciceCheck(not hasToReceiveInitData() && not hasToSendInitData(), "advance()",
	       "initializeData() needs to be called before advance if data has to be initialized!");
  setHasDataBeenExchanged(false);
  setIsCouplingTimestepComplete(false);
  double eps = std::pow(10.0, -1 * getValidDigits());
  
  if (tarch::la::equals(getThisTimestepRemainder(), 0.0, eps)) {
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
      receiveAndSetDt();
      receiveData(getCommunication());
      getCommunication()->finishReceivePackage();
    }
    setHasDataBeenExchanged(true);
    setComputedTimestepPart(0.0);
  }
}

}} // namespace precice, cplscheme
