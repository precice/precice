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
        secondParticipant,localParticipant,communication,maxIterations,dtMethod)
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
    setupDataMatrices(); // Reserve memory and initialize data with zero
    if (_postProcessing.get() != NULL){
      _postProcessing->initialize(getSendData()); // Reserve memory, initialize
    }
  }
  else if (_postProcessing.get() != NULL){ //first participant
    int dataID = _postProcessing->getDataID();
    //TODO
    preciceCheck(getSendData(dataID) == NULL, "initialize()",
                 "A post-processing can be defined for data of second "
                 << "participant only!");
  }

  requireAction(constants::actionWriteIterationCheckpoint());

  // Determine data initialization TODO kommt wahrscheinlich raus
  bool doesReceiveData = not _doesFirstStep;

  // If the second participant initializes data, the first receive for the
  // second participant is done in initializeData() instead of initialize().
  foreach (DataMap::value_type & pair, getSendData()){
    if (pair.second.initialize){
      preciceCheck(not _doesFirstStep, "initialize()",
                   "Only second participant can initialize data!");
      requireAction(constants::actionWriteInitialData());
      preciceDebug("Initialized data to be written");
      doesReceiveData = false;
      break;
    }
  }
  // If the second participant initializes data, the first receive for the first
  // participant is done in initialize() instead of advance().
  // TODO kommt wohl so ähnlich ins initializeData, hier nicht!
  foreach (DataMap::value_type & pair, getReceiveData()){
    if (pair.second.initialize){
      preciceCheck(_doesFirstStep, "initialize()",
                   "Only first participant can receive initial data!");
      preciceDebug("Initialized data to be received");
      doesReceiveData = true;
    }
  }

  if (doesReceiveData && isCouplingOngoing()){
    preciceDebug("Receiving data");
    _communication->startReceivePackage(0);
    if (_participantReceivesDt){ //geht im parallelen sowieso nicht, variable ins
      // neue implicit runter ziehen
      double dt = UNDEFINED_TIMESTEP_LENGTH;
      _communication->receive(dt, 0);
      preciceDebug("received timestep length of " << dt);
      assertion(not tarch::la::equals(dt, UNDEFINED_TIMESTEP_LENGTH));
      setTimestepLength(dt);
      //setMaxLengthNextTimestep(dt);
    }
    receiveData(_communication);
    _communication->finishReceivePackage();
    setHasDataBeenExchanged(true);
  } // TODO bis hier wahrscheinlich raus


  initializeTXTWriters();
  setIsInitialized(true);
}

void ParallelImplicitCouplingScheme:: initializeData()
{

  preciceTrace("initializeData()");
  preciceCheck(isInitialized(), "initializeData()",
               "initializeData() can be called after initialize() only!");
  preciceCheck(isActionRequired(constants::actionWriteInitialData()),
               "initializeData()", "Not required data initialization!");

  //TODO wird nur für 2. gemacht
  foreach (DataMap::value_type & pair, getSendData()){
    utils::DynVector& oldValues = pair.second.oldValues.column(0);
    oldValues = *pair.second.values;

    // For extrapolation, treat the initial value as old timestep value
    pair.second.oldValues.shiftSetFirst(*pair.second.values);
  }

  //TODO gleiche Schleife über receive daten, da die auch gespeichert werden müssen


  //TODO checken ob initialized werden soll, 1. und 2.
  //TODO F: send, receive, S: receive, send
  // The second participant sends the initialized data to the first particpant
  // here, which receives the data on call of initialize().
  sendData(_communication); // TODO sendet immer alle Daten
  _communication->startReceivePackage(0);
  // This receive replaces the receive in initialize().
  receiveData(_communication);
  _communication->finishReceivePackage();
  setHasDataBeenExchanged(true);

  performedAction(constants::actionWriteInitialData());
}


void ParallelImplicitCouplingScheme:: advance()
{
  preciceTrace2("advance()", getTimesteps(), getTime());
  checkCompletenessRequiredActions();
  setHasDataBeenExchanged(false);
  setIsCouplingTimestepComplete(false);
  double eps = std::pow(10.0, -1 * getValidDigits());
  bool convergence = false;
  if (tarch::la::equals(getThisTimestepRemainder(), 0.0, eps)){
    preciceDebug("Computed full length of iteration");
    if (_doesFirstStep){ //FIrst participant
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

      convergence = measureConvergence(); //TODO nochmal checken ob nicht gleicher fall
      //wie bei extrapolate

      assertion2((getSubIteration() <= _maxIterations) || (_maxIterations == -1),
                 getSubIteration(), _maxIterations);
      // Stop, when maximal iteration count (given in config) is reached
      if (getSubIteration() == _maxIterations-1){
        convergence = true;
      }
      if (convergence){
        if (_postProcessing.get() != NULL){
          //TODO embraceData in dieser Klasse (receive und send in eine Map mergen)
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
          //TODO muss auch über alle daten gehen -> extrapolateData so umschreiben
          // im oberen, dass argument hat
          extrapolateData(); // Also stores data
        }
        else { // Store data for conv. measurement, post-processing, or extrapolation
          foreach (DataMap::value_type& pair, getSendData()){
            if (pair.second.oldValues.size() > 0){
              pair.second.oldValues.column(0) = *pair.second.values;
            }
          }
          foreach (DataMap::value_type& pair, getReceiveData()){
            if (pair.second.oldValues.size() > 0){
              pair.second.oldValues.column(0) = *pair.second.values;
            }
          }
        }
        sendData(_communication);
        _communication->finishSendPackage();

      }
      else {
        _communication->finishSendPackage();
      }
    }

    // both participants
    if (not convergence){
      preciceDebug("No convergence achieved");
      requireAction(constants::actionReadIterationCheckpoint());
      setSubIteration(getSubIteration() + 1);
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
      setSubIteration(0);
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
    _iterationToPlot = getSubIteration();
  }
  else {
    _iterationToPlot++;
  }
}


}} // namespace precice, cplscheme
