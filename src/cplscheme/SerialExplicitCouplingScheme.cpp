// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "SerialExplicitCouplingScheme.hpp"
#include "Constants.hpp"
#include "mesh/SharedPointer.hpp"
#include "com/Communication.hpp"
#include "com/SharedPointer.hpp"
#include "tarch/plotter/globaldata/TXTTableWriter.h"

#include "impl/PostProcessing.hpp"

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
			 secondParticipant,localParticipant,communication,dtMethod),
  _postProcessing(),
  _extrapolationOrder(0),
  _convergenceMeasures()  
{}

// SerialExplicitCouplingScheme::initialize and SerialImplicitCouplingScheme::initialize
// are identical now
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
  // This currently does not fail, though description suggests it should in some cases for explicit coupling. 
  preciceCheck(not getSendData().empty(), "initialize()",
               "No send data configured! Use explicit scheme for one-way coupling.");
  setTime(startTime);
  setTimesteps(startTimestep);

  if (not doesFirstStep()){
    if (not _convergenceMeasures.empty()) {
      setupConvergenceMeasures(); // needs _couplingData configured
      setupDataMatrices(getSendData()); // Reserve memory and initialize data with zero
    }
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

  // This test is valid, if only implicit schemes have convergence measures.
  // It currently holds, we will maybe find something better
  if (not _convergenceMeasures.empty()) {
      requireAction(constants::actionWriteIterationCheckpoint());
  }
  
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

void SerialExplicitCouplingScheme::setupDataMatrices(DataMap& data)
{
  preciceTrace("setupDataMatrices()");
  preciceDebug("Data size: " << data.size());
  // Reserve storage for convergence measurement of send and receive data values
  foreach (ConvergenceMeasure& convMeasure, _convergenceMeasures){
    assertion(convMeasure.data != NULL);
    if (convMeasure.data->oldValues.cols() < 1){
      convMeasure.data->oldValues.append(CouplingData::DataMatrix(
          convMeasure.data->values->size(), 1, 0.0));
    }
  }
  // Reserve storage for extrapolation of data values
  if (_extrapolationOrder > 0){
    foreach (DataMap::value_type& pair, data){
      int cols = pair.second->oldValues.cols();
      preciceDebug("Add cols: " << pair.first << ", cols: " << cols);
      assertion1(cols <= 1, cols);
      pair.second->oldValues.append(CouplingData::DataMatrix(
          pair.second->values->size(), _extrapolationOrder + 1 - cols, 0.0));
    }
  }
}

void SerialExplicitCouplingScheme::setupConvergenceMeasures()
{
  preciceTrace("setupConvergenceMeasures()");
  assertion(not doesFirstStep());
  preciceCheck(not _convergenceMeasures.empty(), "setupConvergenceMeasures()",
      "At least one convergence measure has to be defined for "
      << "an implicit coupling scheme!");
  foreach (ConvergenceMeasure& convMeasure, _convergenceMeasures){
    int dataID = convMeasure.dataID;
    if ((getSendData(dataID) != NULL)){
      convMeasure.data = getSendData(dataID);
    }
    else {
      convMeasure.data = getReceiveData(dataID);
      assertion(convMeasure.data != NULL);
    }
  }
}

  
void SerialExplicitCouplingScheme::initializeTXTWriters()
{
  // NoOp
}



}} // namespace precice, cplscheme
