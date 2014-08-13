#include "ParallelCouplingScheme.hpp"
#include "impl/PostProcessing.hpp"
#include "com/Communication.hpp"
#include "com/SharedPointer.hpp"

namespace precice {
namespace cplscheme {

tarch::logging::Log ParallelCouplingScheme::_log("precice::cplscheme::ParallelCouplingScheme" );

ParallelCouplingScheme::ParallelCouplingScheme
(
  double                maxTime,
  int                   maxTimesteps,
  double                timestepLength,
  int                   validDigits,
  const std::string&    firstParticipant,
  const std::string&    secondParticipant,
  const std::string&    localParticipant,
  com::PtrCommunication communication,
  constants::TimesteppingMethod dtMethod,
  CouplingMode          cplMode,
  int                   maxIterations)
  :
  BaseCouplingScheme(maxTime,maxTimesteps,timestepLength,validDigits,firstParticipant,
                     secondParticipant,localParticipant,communication,maxIterations,dtMethod),
  _allData ()
{
  _couplingMode = cplMode;
  // Coupling mode must be either Explicit or Implicit when using SerialCouplingScheme.
  assertion(_couplingMode != Undefined);
  if (_couplingMode == Explicit) {
    assertion(maxIterations == 1);
  }
}

void ParallelCouplingScheme::initialize
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
  if (_couplingMode == Implicit) {
    if (not doesFirstStep()) { // second participant
      setupConvergenceMeasures(); // needs _couplingData configured
      mergeData(); // merge send and receive data for all pp calls
      setupDataMatrices(getAllData()); // Reserve memory and initialize data with zero
      if (getPostProcessing().get() != NULL) {
        preciceCheck(getPostProcessing()->getDataIDs().size()==2 ,"initialize()",
                     "For parallel coupling, the number of coupling data vectors has to be 2, not: "
                     << getPostProcessing()->getDataIDs().size());
        getPostProcessing()->initialize(getAllData()); // Reserve memory, initialize
      }
    }

    requireAction(constants::actionWriteIterationCheckpoint());
    initializeTXTWriters();
  }

  foreach (DataMap::value_type & pair, getSendData()) {
    if (pair.second->initialize) {
      setHasToSendInitData(true);
      break;
    }
  }
  foreach (DataMap::value_type & pair, getReceiveData()) {
    if (pair.second->initialize) {
      setHasToReceiveInitData(true);
      break;
    }
  }

  if (hasToSendInitData()) {
    requireAction(constants::actionWriteInitialData());
  }

  setIsInitialized(true);
}

void ParallelCouplingScheme::initializeData()
{
  preciceTrace("initializeData()");
  preciceCheck(isInitialized(), "initializeData()",
               "initializeData() can be called after initialize() only!");

  if (not hasToSendInitData() && not hasToReceiveInitData()) {
    preciceInfo("initializeData()", "initializeData is skipped since no data has to be initialized");
    return;
  }

  preciceCheck(not (hasToSendInitData() && isActionRequired(constants::actionWriteInitialData())),
               "initializeData()", "InitialData has to be written to preCICE before calling initializeData()");

  setHasDataBeenExchanged(false);

  // F: send, receive, S: receive, send
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

      // second participant has to save values for extrapolation
      if (_couplingMode == Implicit and getExtrapolationOrder() > 0){
        foreach (DataMap::value_type & pair, getReceiveData()){
          utils::DynVector& oldValues = pair.second->oldValues.column(0);
          oldValues = *pair.second->values;
          // For extrapolation, treat the initial value as old timestep value
          pair.second->oldValues.shiftSetFirst(*pair.second->values);
        }
      }
    }
    if (hasToSendInitData()) {
      if (_couplingMode == Implicit and getExtrapolationOrder() > 0) {
        foreach (DataMap::value_type & pair, getSendData()) {
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

  // in order to check in advance if initializeData has been called (if necessary)
  setHasToSendInitData(false);
  setHasToReceiveInitData(false);
}

void ParallelCouplingScheme::advance()
{
  if (_couplingMode == Explicit) {
    explicitAdvance();
  }
  else if (_couplingMode == Implicit) {
    implicitAdvance();
  }
}

void ParallelCouplingScheme::explicitAdvance()
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
      receiveAndSetDt();
      receiveData(getCommunication());
      getCommunication()->finishReceivePackage();
      setHasDataBeenExchanged(true);
    }
    else { //second participant
      preciceDebug("Receiving data...");
      getCommunication()->startReceivePackage(0);
      receiveAndSetDt();
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

void ParallelCouplingScheme::implicitAdvance()
{
  preciceTrace2("advance()", getTimesteps(), getTime());
  checkCompletenessRequiredActions();

  preciceCheck(!hasToReceiveInitData() && !hasToSendInitData(), "advance()",
               "initializeData() needs to be called before advance if data has to be initialized!");

  setHasDataBeenExchanged(false);
  setIsCouplingTimestepComplete(false);
  double eps = std::pow(10.0, -1 * getValidDigits());
  bool convergence = false;
  if (tarch::la::equals(getThisTimestepRemainder(), 0.0, eps)) {
    // bis hier identisch mit Explicit, von vorne
    preciceDebug("Computed full length of iteration");
    if (doesFirstStep()) { //First participant
      getCommunication()->startSendPackage(0);
      sendData(getCommunication());
      getCommunication()->finishSendPackage();
      getCommunication()->startReceivePackage(0);
      getCommunication()->receive(convergence, 0);
      if (convergence) {
        timestepCompleted();
      }
      if (isCouplingOngoing()) {
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
      if (getIterations() == getMaxIterations()-1) {
        convergence = true;
      }
      if (convergence) {
        if (getPostProcessing().get() != NULL) {
          getPostProcessing()->iterationsConverged(getAllData());
        }
        newConvergenceMeasurements();
        timestepCompleted();
      }
      else if (getPostProcessing().get() != NULL) {
        getPostProcessing()->performPostProcessing(getAllData());
      }
      getCommunication()->startSendPackage(0);
      getCommunication()->send(convergence, 0);

      if (isCouplingOngoing()) {
        if (convergence && (getExtrapolationOrder() > 0)){
          extrapolateData(getAllData()); // Also stores data
        }
        else { // Store data for conv. measurement, post-processing, or extrapolation
          foreach (DataMap::value_type& pair, getSendData()) {
            if (pair.second->oldValues.size() > 0){
              pair.second->oldValues.column(0) = *pair.second->values;
            }
          }
          foreach (DataMap::value_type& pair, getReceiveData()) {
            if (pair.second->oldValues.size() > 0) {
              pair.second->oldValues.column(0) = *pair.second->values;
            }
          }
        }
        sendData(getCommunication());
      }
      getCommunication()->finishSendPackage();
    }

    // both participants
    if (not convergence) {
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
  if (not convergence) {
    setTimestepToPlot(getTimesteps());
    setTimeToPlot(getTime());
    setIterationToPlot(getIterations());
  }
  else {
    increaseIterationToPlot();
  }
}



void ParallelCouplingScheme::mergeData()
{
  preciceTrace("mergeData()");
  assertion1(!doesFirstStep(), "Only the second participant should do the post processing." );
  assertion1(_allData.empty(), "This function should only be called once.");
  _allData.insert(getSendData().begin(), getSendData().end());
  _allData.insert(getReceiveData().begin(), getReceiveData().end());
}

}}
