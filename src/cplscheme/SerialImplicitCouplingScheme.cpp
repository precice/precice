// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "SerialImplicitCouplingScheme.hpp"
#include "impl/PostProcessing.hpp"
#include "impl/ConvergenceMeasure.hpp"
#include "Constants.hpp"
#include "mesh/SharedPointer.hpp"
#include "com/Communication.hpp"
#include "com/SharedPointer.hpp"
#include "tarch/plotter/globaldata/TXTTableWriter.h"

namespace precice {
namespace cplscheme {

tarch::logging::Log SerialImplicitCouplingScheme::
    _log("precice::cplscheme::SerialImplicitCouplingScheme" );

SerialImplicitCouplingScheme:: SerialImplicitCouplingScheme
(
  double                maxTime,
  int                   maxTimesteps,
  double                timestepLength,
  int                   validDigits,
  const std::string&    firstParticipant,
  const std::string&    secondParticipant,
  const std::string&    localParticipant,
  com::PtrCommunication communication,
  int                   maxIterations,
  constants::TimesteppingMethod dtMethod )
:
  ImplicitCouplingScheme(maxTime,maxTimesteps,timestepLength,validDigits,firstParticipant,
        secondParticipant,localParticipant,communication,maxIterations,dtMethod)
{}

SerialImplicitCouplingScheme:: ~SerialImplicitCouplingScheme()
{}



void SerialImplicitCouplingScheme:: initialize
(
  double startTime,
  int    startTimestep )
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
  if (not doesFirstStep()){
    setupConvergenceMeasures(); // needs _couplingData configured
    setupDataMatrices(getSendData()); // Reserve memory and initialize data with zero
    if (getPostProcessing().get() != NULL){
      preciceCheck(getPostProcessing()->getDataIDs().size()==1 ,"initialize()",
                    "For serial coupling, the number of coupling data vectors has to be 1");
      getPostProcessing()->initialize(getSendData()); // Reserve memory, initialize
    }
  }
  else if (getPostProcessing().get() != NULL){
    int dataID = *(getPostProcessing()->getDataIDs().begin());
    preciceCheck(getSendData(dataID) == NULL, "initialize()",
                 "In case of serial coupling, post-processing can be defined for "
                 << "data of second participant only!");
  }

  requireAction(constants::actionWriteIterationCheckpoint());




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

  initializeTXTWriters();
  setIsInitialized(true);
  //tmp debug, here still the right values
      foreach (DataMap::value_type & pair, getReceiveData()){
            utils::DynVector& values = *pair.second->values;
            preciceDebug("End initialize, New Values: " << values);
      }
}

// SerialExplicitCouplingScheme::initializeData and SerialImplicitCouplingScheme::initializeData
// are identical now
void SerialImplicitCouplingScheme:: initializeData()
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


void SerialImplicitCouplingScheme:: advance()
{
  preciceTrace2("advance()", getTimesteps(), getTime());
  //tmp debug -> here wrong values
            foreach (DataMap::value_type & pair, getReceiveData()){
                  utils::DynVector& values = *pair.second->values;
                  preciceDebug("Begin advance, New Values: " << values);
            }
  checkCompletenessRequiredActions();

  preciceCheck(not hasToReceiveInitData() && not hasToSendInitData(), "advance()",
     "initializeData() needs to be called before advance if data has to be initialized!");

  setHasDataBeenExchanged(false);
  setIsCouplingTimestepComplete(false);
  double eps = std::pow(10.0, -1 * getValidDigits());
  bool convergence = false;
  if (tarch::la::equals(getThisTimestepRemainder(), 0.0, eps)){
    preciceDebug("Computed full length of iteration");
    if (doesFirstStep()){
      getCommunication()->startSendPackage(0);
      if (participantSetsDt()){
        preciceDebug("sending timestep length of " << getComputedTimestepPart());
        getCommunication()->send(getComputedTimestepPart(), 0);
      }
      sendData(getCommunication());
      getCommunication()->finishSendPackage();
      getCommunication()->startReceivePackage(0);
      getCommunication()->receive(convergence, 0);
      if (convergence){
        timestepCompleted();
      }
      if (isCouplingOngoing()){
        receiveData(getCommunication());
      }
      getCommunication()->finishReceivePackage();
    }
    else {
      convergence = measureConvergence();
      assertion2((getIterations() <= getMaxIterations()) || (getMaxIterations() == -1),
                 getIterations(), getMaxIterations());
      // Stop, when maximal iteration count (given in config) is reached
      if (getIterations() == getMaxIterations()-1){
        convergence = true;
      }
      if (convergence){
        if (getPostProcessing().get() != NULL){
          getPostProcessing()->iterationsConverged(getSendData());
        }
        newConvergenceMeasurements();
        timestepCompleted();
      }
      else if (getPostProcessing().get() != NULL){
        getPostProcessing()->performPostProcessing(getSendData());
      }
      getCommunication()->startSendPackage(0);
      getCommunication()->send(convergence, 0);
      if (isCouplingOngoing()){
        if (convergence && (getExtrapolationOrder() > 0)){
          extrapolateData(getSendData()); // Also stores data
        }
        else { // Store data for conv. measurement, post-processing, or extrapolation
          foreach (DataMap::value_type& pair, getSendData()){
            if (pair.second->oldValues.size() > 0){
              pair.second->oldValues.column(0) = *pair.second->values;
            }
          }
          foreach (DataMap::value_type& pair, getReceiveData()){
            if (pair.second->oldValues.size() > 0){
              pair.second->oldValues.column(0) = *pair.second->values;
            }
          }
        }
        sendData(getCommunication());
        getCommunication()->finishSendPackage();
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
      else {
        getCommunication()->finishSendPackage();
      }
    }

    if (not convergence){
      preciceDebug("No convergence achieved");
      requireAction(constants::actionReadIterationCheckpoint());
      increaseIterations();
      increaseTotalIterations();
      // The computed timestep part equals the timestep length, since the
      // timestep remainder is zero. Subtract the timestep length do another
      // coupling iteration.
      assertion(tarch::la::greater(getComputedTimestepPart(), 0.0));
      setTime(getTime() - getComputedTimestepPart());
    }
    else {
      preciceDebug("Convergence achieved");
      getIterationsWriter().writeData("Timesteps", getTimesteps());
      getIterationsWriter().writeData("Total Iterations", getTotalIterations());
      getIterationsWriter().writeData("Iterations", getIterations());
      int converged = getIterations() < getMaxIterations() ? 1 : 0;
      getIterationsWriter().writeData("Convergence", converged);
      setIterations(0);
    }
    setHasDataBeenExchanged(true);
    setComputedTimestepPart(0.0);
  }//subcycling completed

  // When the iterations of one timestep are converged, the old time, timesteps,
  // and iteration should be plotted, and not the 0th of the new timestep. Thus,
  // the plot values are only updated when no convergence was achieved.
  if (not convergence){
    setTimestepToPlot(getTimesteps());
    setTimeToPlot(getTime());
    setIterationToPlot(getIterations());
  }
  else {
    increaseIterationToPlot();
  }

}

}} // namespace precice, cplscheme
