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
    _log("precice::cplscheme::ImplicitCouplingScheme" );

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
  _allData (),
  _hasToSendInitData (false),
  _hasToReceiveInitData (false)
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
  assertion(_communication->isConnected());
  preciceCheck(not getSendData().empty(), "initialize()",
               "No send data configured!");
  setTime(startTime);
  setTimesteps(startTimestep);
  if (not _doesFirstStep){ // second participant
    setupConvergenceMeasures(); // needs _couplingData configured
    mergeData(); // merge send and receive data for all pp calls
    setupDataMatrices(getAllData()); // Reserve memory and initialize data with zero
    if (_postProcessing.get() != NULL){
      _postProcessing->initialize(getAllData()); // Reserve memory, initialize
    }
  }


  initializeTXTWriters();


  foreach (DataMap::value_type & pair, getSendData()){
    if (pair.second->initialize){
      _hasToSendInitData = true;
      break;
    }
  }
  foreach (DataMap::value_type & pair, getReceiveData()){
    if (pair.second->initialize){
      _hasToReceiveInitData = true;
      break;
    }
  }

  if(_hasToSendInitData){
    requireAction(constants::actionWriteInitialData());
  }

  setIsInitialized(true);

}

void ParallelImplicitCouplingScheme:: initializeData()
{
  preciceTrace("initializeData()");
  preciceCheck(isInitialized(), "initializeData()",
     "initializeData() can be called after initialize() only!");

  if(not _hasToSendInitData && not _hasToReceiveInitData){
    preciceInfo("initializeData()", "initializeData is skipped since no data has to be initialized");
    return;
  }

  preciceCheck(not (_hasToSendInitData && isActionRequired(constants::actionWriteInitialData())),
     "initializeData()", "InitialData has to be written to preCICE before calling initializeData()");

  setHasDataBeenExchanged(false);


  //F: send, receive, S: receive, send
  if(_doesFirstStep){
    if(_hasToSendInitData){
      _communication->startSendPackage(0);
      sendData(_communication);
      _communication->finishSendPackage();
    }
    if(_hasToReceiveInitData){
      _communication->startReceivePackage(0);
      receiveData(_communication);
      _communication->finishReceivePackage();

      setHasDataBeenExchanged(true);
    }

  }

  else{ // second participant
    if(_hasToReceiveInitData){


      _communication->startReceivePackage(0);
      receiveData(_communication);
      _communication->finishReceivePackage();

      setHasDataBeenExchanged(true);


      // second participant has to save values for extrapolation
      foreach (DataMap::value_type & pair, getReceiveData()){
        utils::DynVector& oldValues = pair.second->oldValues.column(0);
        oldValues = *pair.second->values;
        // For extrapolation, treat the initial value as old timestep value
        pair.second->oldValues.shiftSetFirst(*pair.second->values);
        preciceDebug("Shift columns for receive data " << pair.first);
        preciceDebug("Columns for receive data " << pair.second->oldValues.cols());
      }

    }
    if(_hasToSendInitData){

      foreach (DataMap::value_type & pair, getSendData()){
        utils::DynVector& oldValues = pair.second->oldValues.column(0);
        oldValues = *pair.second->values;
        // For extrapolation, treat the initial value as old timestep value
        pair.second->oldValues.shiftSetFirst(*pair.second->values);
        preciceDebug("old Values in initializeData: " << oldValues);
        preciceDebug("Shift columns for send data " << pair.first);
        preciceDebug("columns for send data " << pair.second->oldValues.cols());
      }

      _communication->startSendPackage(0);
      sendData(_communication);
      _communication->finishSendPackage();


    }

  }

  //in order to check in advance if initializeData has been called (if necessary)
  _hasToSendInitData = false;
  _hasToReceiveInitData = false;

}


void ParallelImplicitCouplingScheme:: advance()
{
  preciceTrace2("advance()", getTimesteps(), getTime());
  checkCompletenessRequiredActions();

  preciceCheck(!_hasToReceiveInitData && !_hasToSendInitData, "advance()",
     "initializeData() needs to be called before advance if data has to be initialized!");

  setHasDataBeenExchanged(false);
  setIsCouplingTimestepComplete(false);
  double eps = std::pow(10.0, -1 * getValidDigits());
  bool convergence = false;
  if (tarch::la::equals(getThisTimestepRemainder(), 0.0, eps)){
    preciceDebug("Computed full length of iteration");
    if (_doesFirstStep){ //First participant
      _communication->startSendPackage(0);
      sendData(_communication);
      _communication->finishSendPackage();
      _communication->startReceivePackage(0);
      _communication->receive(convergence, 0);
      if (convergence){
        timestepCompleted();
      }
      if (isCouplingOngoing()){
        receiveData(_communication);
      }
      _communication->finishReceivePackage();
    }
    else { // second participant

      _communication->startReceivePackage(0);
      receiveData(_communication);
      _communication->finishReceivePackage();

      convergence = measureConvergence();

      //TODO Debug
      if(not _doesFirstStep){
        foreach (DataMap::value_type& pair, getSendData()){
          preciceDebug("after MC: " << pair.second->oldValues.column(0));
        }
      }
      if(not _doesFirstStep){
        foreach (DataMap::value_type& pair, getAllData()){
          preciceDebug("after MC (ALL DATA): " << pair.second->oldValues.column(0));
        }
      }

      assertion2((_iterations <= _maxIterations) || (_maxIterations == -1),
                    _iterations, _maxIterations);
      // Stop, when maximal iteration count (given in config) is reached
      if (_iterations == _maxIterations-1){
        convergence = true;
      }
      if (convergence){
        if (_postProcessing.get() != NULL){
          _postProcessing->iterationsConverged(getAllData());
        }
        newConvergenceMeasurements();
        timestepCompleted();
      }
      else if (_postProcessing.get() != NULL){
        _postProcessing->performPostProcessing(getAllData());
      }
      _communication->startSendPackage(0);
      _communication->send(convergence, 0);

      if (isCouplingOngoing()){
        if (convergence && (_extrapolationOrder > 0)){
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
        sendData(_communication);
      }

      _communication->finishSendPackage();

    }

    // both participants
    if (not convergence){
      preciceDebug("No convergence achieved");
      requireAction(constants::actionReadIterationCheckpoint());
      _iterations++;
      _totalIterations++;
      // The computed timestep part equals the timestep length, since the
      // timestep remainder is zero. Subtract the timestep length do another
      // coupling iteration.
      assertion(tarch::la::greater(getComputedTimestepPart(), 0.0));
      setTime(getTime() - getComputedTimestepPart());
    }
    else {
      preciceDebug("Convergence achieved");
      _iterationsWriter.writeData("Timesteps", getTimesteps());
      _iterationsWriter.writeData("Total Iterations", _totalIterations);
      _iterationsWriter.writeData("Iterations", getSubIteration());
      int converged = getSubIteration() < _maxIterations ? 1 : 0;
      _iterationsWriter.writeData("Convergence", converged);
      _iterations = 0;
    }
    setHasDataBeenExchanged(true);
    setComputedTimestepPart(0.0);
  } // subcycling complete

  // When the iterations of one timestep are converged, the old time, timesteps,
  // and iteration should be plotted, and not the 0th of the new timestep. Thus,
  // the plot values are only updated when no convergence was achieved.
  if (not convergence){
    _timestepToPlot = getTimesteps();
    _timeToPlot = getTime();
    _iterationToPlot = _iterations;
  }
  else {
    _iterationToPlot++;
  }
}

void ParallelImplicitCouplingScheme:: mergeData()
{
  preciceTrace("mergeData()");
  assertion1(!_doesFirstStep, "Only the second participant should do the pp." );
  assertion1(_allData.empty(), "This function should only be called once.");
  _allData.insert(getSendData().begin(),getSendData().end());
  _allData.insert(getReceiveData().begin(),getReceiveData().end());

}

}} // namespace precice, cplscheme
