// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "ParallelExplicitCouplingScheme.hpp"
#include "impl/PostProcessing.hpp"
#include "impl/ConvergenceMeasure.hpp"
#include "Constants.hpp"
#include "com/Communication.hpp"
#include "com/SharedPointer.hpp"
#include "io/TXTWriter.hpp"
#include "io/TXTReader.hpp"
#include "tarch/plotter/globaldata/TXTTableWriter.h"

namespace precice {
namespace cplscheme {

tarch::logging::Log ParallelExplicitCouplingScheme::
_log("precice::cplscheme::ParallelExplicitCouplingScheme" );

ParallelExplicitCouplingScheme:: ParallelExplicitCouplingScheme
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
  ParallelCouplingScheme(maxTime,maxTimesteps,timestepLength,validDigits,firstParticipant,
			 secondParticipant,localParticipant,communication, 1, dtMethod)
{
  couplingMode = Explicit; 
}

void ParallelExplicitCouplingScheme:: initializeData()
{
  preciceTrace("initializeData()");
  preciceCheck(isInitialized(), "initializeData()",
	       "initializeData() can be called after initialize() only!");

  if(not hasToSendInitData() && not hasToReceiveInitData()){
    preciceInfo("initializeData()", "initializeData is skipped since no data has to be initialized");
    return;
  }

  preciceCheck(not (hasToSendInitData() && isActionRequired(constants::actionWriteInitialData())),
	       "initializeData()", "InitialData has to be written to preCICE before calling initializeData()");

  setHasDataBeenExchanged(false);

  //F: send, receive, S: receive, send
  if (doesFirstStep()) {
    if (hasToSendInitData()) {
      getCommunication()->startSendPackage(0);
      sendData(getCommunication());
      getCommunication()->finishSendPackage();
    }
    if (hasToReceiveInitData()) {
      getCommunication()->startReceivePackage(0);
      receiveData(getCommunication());
      getCommunication()->finishReceivePackage();
      setHasDataBeenExchanged(true);
    }
  }

  else { // second participant
    if (hasToReceiveInitData()) {
      getCommunication()->startReceivePackage(0);
      receiveData(getCommunication());
      getCommunication()->finishReceivePackage();
      setHasDataBeenExchanged(true);
    }
    if (hasToSendInitData()) {
      getCommunication()->startSendPackage(0);
      sendData(getCommunication());
      getCommunication()->finishSendPackage();
    }
  }

  //in order to check in advance if initializeData has been called (if necessary)
  setHasToSendInitData(false);
  setHasToReceiveInitData(false);
}

void ParallelExplicitCouplingScheme:: advance()
{
  preciceTrace("advance()");
  checkCompletenessRequiredActions();
  preciceCheck(!hasToReceiveInitData() && !hasToSendInitData(), "advance()",
	       "initializeData() needs to be called before advance if data has to be initialized!");
  setHasDataBeenExchanged(false);
  setIsCouplingTimestepComplete(false);

  double eps = std::pow(10.0, -1 * getValidDigits());
  if (tarch::la::equals(getThisTimestepRemainder(), 0.0, eps)) {
    setIsCouplingTimestepComplete(true);
    setTimesteps(getTimesteps() + 1);

    if (doesFirstStep()) {
      preciceDebug("Sending data...");
      getCommunication()->startSendPackage(0);
      if (participantSetsDt()) {
        getCommunication()->send(getComputedTimestepPart(), 0);
      }
      sendData(getCommunication());
      getCommunication()->finishSendPackage();

      preciceDebug("Receiving data...");
      getCommunication()->startReceivePackage(0);
      if (participantReceivesDt()) {
        double dt = UNDEFINED_TIMESTEP_LENGTH;
        getCommunication()->receive(dt, 0);
        assertion(not tarch::la::equals(dt, UNDEFINED_TIMESTEP_LENGTH));
        setTimestepLength(dt);
      }
      receiveData(getCommunication());
      getCommunication()->finishReceivePackage();
      setHasDataBeenExchanged(true);
    }

    else{ //second participant
      preciceDebug("Receiving data...");
      getCommunication()->startReceivePackage(0);
      if (participantReceivesDt()) {
        double dt = UNDEFINED_TIMESTEP_LENGTH;
        getCommunication()->receive(dt, 0);
        assertion(not tarch::la::equals(dt, UNDEFINED_TIMESTEP_LENGTH));
        setTimestepLength(dt);
      }
      receiveData(getCommunication());
      getCommunication()->finishReceivePackage();
      setHasDataBeenExchanged(true);

      preciceDebug("Sending data...");
      getCommunication()->startSendPackage(0);
      if (participantSetsDt()) {
        getCommunication()->send(getComputedTimestepPart(), 0);
      }
      sendData(getCommunication());
      getCommunication()->finishSendPackage();
    }

    //both participants
    setComputedTimestepPart(0.0);
  }
}

}} // namespace precice, cplscheme
