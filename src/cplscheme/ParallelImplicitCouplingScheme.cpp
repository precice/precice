// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "ParallelImplicitCouplingScheme.hpp"
#include "impl/PostProcessing.hpp"
#include "impl/ConvergenceMeasure.hpp"
#include "Constants.hpp"
#include "mesh/SharedPointer.hpp"
#include "com/Communication.hpp"
#include "com/SharedPointer.hpp"
#include "io/TXTWriter.hpp"
#include "io/TXTReader.hpp"
#include "tarch/plotter/globaldata/TXTTableWriter.h"
#include <limits>

namespace precice {
namespace cplscheme {

tarch::logging::Log ParallelImplicitCouplingScheme::
    _log("precice::cplscheme::ParallelImplicitCouplingScheme" );

ParallelImplicitCouplingScheme:: ParallelImplicitCouplingScheme
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
        secondParticipant,localParticipant,communication,maxIterations,dtMethod),
  _allData ()
{}

ParallelImplicitCouplingScheme:: ~ParallelImplicitCouplingScheme()
{}



void ParallelImplicitCouplingScheme:: initialize
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
               "No send data configured!");
  setTime(startTime);
  setTimesteps(startTimestep);
  if (not doesFirstStep()){ // second participant
    setupConvergenceMeasures(); // needs _couplingData configured
    mergeData(); // merge send and receive data for all pp calls
    setupDataMatrices(getAllData()); // Reserve memory and initialize data with zero
    if (getPostProcessing().get() != NULL){
      preciceCheck(getPostProcessing()->getDataIDs().size()==2 ,"initialize()",
              "For parallel coupling, the number of coupling data vectors has to be 2, not: "
              << getPostProcessing()->getDataIDs().size());
      getPostProcessing()->initialize(getAllData()); // Reserve memory, initialize
    }
  }

  requireAction(constants::actionWriteIterationCheckpoint());
  initializeTXTWriters();

  foreach (DataMap::value_type & pair, getSendData()){
    if (pair.second->initialize){
      setHasToSendInitData(true);
      break;
    }
  }
  foreach (DataMap::value_type & pair, getReceiveData()){
    if (pair.second->initialize){
      setHasToReceiveInitData(true);
      break;
    }
  }

  if(hasToSendInitData()){
    requireAction(constants::actionWriteInitialData());
  }

  setIsInitialized(true);
}

void ParallelImplicitCouplingScheme:: initializeData()
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
  if(doesFirstStep()){
    if(hasToSendInitData()){
      getCommunication()->startSendPackage(0);
      sendData(getCommunication());
      getCommunication()->finishSendPackage();
    }
    if(hasToReceiveInitData()){
      getCommunication()->startReceivePackage(0);
      receiveData(getCommunication());
      getCommunication()->finishReceivePackage();
      setHasDataBeenExchanged(true);
    }
  }

  else{ // second participant
    if(hasToReceiveInitData()){
      getCommunication()->startReceivePackage(0);
      receiveData(getCommunication());
      getCommunication()->finishReceivePackage();
      setHasDataBeenExchanged(true);

      // second participant has to save values for extrapolation
      if (getExtrapolationOrder() > 0){
        foreach (DataMap::value_type & pair, getReceiveData()){
          utils::DynVector& oldValues = pair.second->oldValues.column(0);
          oldValues = *pair.second->values;
          // For extrapolation, treat the initial value as old timestep value
          pair.second->oldValues.shiftSetFirst(*pair.second->values);
        }
      }
    }
    if(hasToSendInitData()){
      if (getExtrapolationOrder() > 0){
        foreach (DataMap::value_type & pair, getSendData()){
          utils::DynVector& oldValues = pair.second->oldValues.column(0);
          oldValues = *pair.second->values;
          // For extrapolation, treat the initial value as old timestep value
          pair.second->oldValues.shiftSetFirst(*pair.second->values);
        }
      }
      getCommunication()->startSendPackage(0);
      sendData(getCommunication());
      getCommunication()->finishSendPackage();
    }
  }

  //in order to check in advance if initializeData has been called (if necessary)
  setHasToSendInitData(false);
  setHasToReceiveInitData(false);
}


void ParallelImplicitCouplingScheme:: advance()
{
  preciceTrace2("advance()", getTimesteps(), getTime());
  checkCompletenessRequiredActions();

  preciceCheck(!hasToReceiveInitData() && !hasToSendInitData(), "advance()",
     "initializeData() needs to be called before advance if data has to be initialized!");

  setHasDataBeenExchanged(false);
  setIsCouplingTimestepComplete(false);
  double eps = std::pow(10.0, -1 * getValidDigits());
  bool convergence = false;
  if (tarch::la::equals(getThisTimestepRemainder(), 0.0, eps)){
    preciceDebug("Computed full length of iteration");
    if (doesFirstStep()){ //First participant
      getCommunication()->startSendPackage(0);
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
    else { // second participant

      getCommunication()->startReceivePackage(0);
      receiveData(getCommunication());
      getCommunication()->finishReceivePackage();

      convergence = measureConvergence();

      assertion2((getIterations() <= getMaxIterations()) || (getMaxIterations() == -1),
                    getIterations(), getMaxIterations());
      // Stop, when maximal iteration count (given in config) is reached
      if (getIterations() == getMaxIterations()-1){
        convergence = true;
      }
      if (convergence){
        if (getPostProcessing().get() != NULL){
          getPostProcessing()->iterationsConverged(getAllData());
        }
        newConvergenceMeasurements();
        timestepCompleted();
      }
      else if (getPostProcessing().get() != NULL){
        getPostProcessing()->performPostProcessing(getAllData());
      }
      getCommunication()->startSendPackage(0);
      getCommunication()->send(convergence, 0);

      if (isCouplingOngoing()){
        if (convergence && (getExtrapolationOrder() > 0)){
          extrapolateData(getAllData()); // Also stores data
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
      }

      getCommunication()->finishSendPackage();

    }

    // both participants
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
  } // subcycling complete

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

void ParallelImplicitCouplingScheme:: mergeData()
{
  preciceTrace("mergeData()");
  assertion1(!doesFirstStep(), "Only the second participant should do the pp." );
  assertion1(_allData.empty(), "This function should only be called once.");
  _allData.insert(getSendData().begin(),getSendData().end());
  _allData.insert(getReceiveData().begin(),getReceiveData().end());

}

}} // namespace precice, cplscheme
