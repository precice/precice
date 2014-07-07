// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "SerialExplicitCouplingScheme.hpp"
#include "Constants.hpp"
#include "mesh/SharedPointer.hpp"
#include "com/Communication.hpp"
#include "com/SharedPointer.hpp"
#include "tarch/plotter/globaldata/TXTTableWriter.h"

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

SerialExplicitCouplingScheme:: ~SerialExplicitCouplingScheme()
{}

void SerialExplicitCouplingScheme:: initialize
(
  double startTime,
  int    startTimestep)
{
  preciceTrace2("initialize()", startTime, startTimestep);
  assertion(not isInitialized());
  assertion1(tarch::la::greaterEquals(startTime, 0.0), startTime);
  assertion1(startTimestep >= 0, startTimestep);
  assertion(getCommunication()->isConnected());
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

// SerialExplicitCouplingScheme::initializeData and SerialImplicitCouplingScheme::initializeData
// are identical now
void SerialExplicitCouplingScheme:: initializeData()
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
    foreach (DataMap::value_type & pair, getSendData()){
      if (pair.second->oldValues.cols() == 0)
	break;
      utils::DynVector& oldValues = pair.second->oldValues.column(0);
      oldValues = *pair.second->values;

      // For extrapolation, treat the initial value as old timestep value
      pair.second->oldValues.shiftSetFirst(*pair.second->values);
    }

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

void SerialExplicitCouplingScheme:: advance()
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




}} // namespace precice, cplscheme
