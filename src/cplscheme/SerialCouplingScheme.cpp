#include "SerialCouplingScheme.hpp"
#include "impl/PostProcessing.hpp"

namespace precice {
namespace cplscheme {

tarch::logging::Log SerialCouplingScheme::_log("precice::cplscheme::SerialCouplingScheme" );

SerialCouplingScheme::SerialCouplingScheme
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
  BaseCouplingScheme(maxTime, maxTimesteps, timestepLength, validDigits, firstParticipant,
                     secondParticipant, localParticipant, communication, maxIterations, dtMethod)
{
  _couplingMode = cplMode;
  // Coupling mode must be either Explicit or Implicit when using SerialCouplingScheme.
  assertion(_couplingMode != Undefined);
  if (_couplingMode == Explicit) {
    assertion(maxIterations == 1);
  }
}

void SerialCouplingScheme::initialize
(
  double startTime,
  int    startTimestep)
{
  preciceTrace2("initialize()", startTime, startTimestep);
  assertion(not isInitialized());
  assertion1(tarch::la::greaterEquals(startTime, 0.0), startTime);
  assertion1(startTimestep >= 0, startTimestep);
  setTime(startTime);
  setTimesteps(startTimestep);
  
  if (_couplingMode == Implicit) {
    preciceCheck(not getSendData().empty(), "initialize()", "No send data configured! Use explicit scheme for one-way coupling.");
    if (not doesFirstStep()) {
      if (not _convergenceMeasures.empty()) {
        setupConvergenceMeasures(); // needs _couplingData configured
        setupDataMatrices(getSendData()); // Reserve memory and initialize data with zero
      }
      if (getPostProcessing().get() != NULL) {
        preciceCheck(getPostProcessing()->getDataIDs().size()==1 ,"initialize()",
                     "For serial coupling, the number of coupling data vectors has to be 1");
        getPostProcessing()->initialize(getSendData()); // Reserve memory, initialize
      }
    }
    else if (getPostProcessing().get() != NULL) {
      int dataID = *(getPostProcessing()->getDataIDs().begin());
      preciceCheck(getSendData(dataID) == NULL, "initialize()",
                   "In case of serial coupling, post-processing can be defined for "
                   << "data of second participant only!");
    }
    requireAction(constants::actionWriteIterationCheckpoint());
  }
    
  foreach (DataMap::value_type & pair, getSendData()){
    if (pair.second->initialize) {
      preciceCheck(not doesFirstStep(), "initialize()",
                   "Only second participant can initialize data!");
      preciceDebug("Initialized data to be written");
      setHasToSendInitData(true);
      break;
    }
  }

  foreach (DataMap::value_type & pair, getReceiveData()){
    if (pair.second->initialize) {
      preciceCheck(doesFirstStep(), "initialize()",
                   "Only first participant can receive initial data!");
      preciceDebug("Initialized data to be received");
      setHasToReceiveInitData(true);
    }
  }
  
  // If the second participant initializes data, the first receive for the
  // second participant is done in initializeData() instead of initialize().
  if (not doesFirstStep() && not hasToSendInitData() && isCouplingOngoing()) {
    preciceDebug("Receiving data");
    startReceivePackage();
    receiveAndSetDt();
    receiveData(getCommunication());
    finishReceivePackage();
    setHasDataBeenExchanged(true);
  }

  if (hasToSendInitData()) {
    requireAction(constants::actionWriteInitialData());
  }
  
  initializeTXTWriters();
  setIsInitialized(true);
}


void SerialCouplingScheme::initializeData()
{
  preciceTrace("initializeData()");
  preciceCheck(isInitialized(), "initializeData()",
               "initializeData() can be called after initialize() only!");

  if (not hasToSendInitData() && not hasToReceiveInitData()) {
    preciceInfo("initializeData()", "initializeData is skipped since no data has to be initialized");
    return;
  }

  preciceDebug("Initializing Data ...");
  
  preciceCheck(not (hasToSendInitData() && isActionRequired(constants::actionWriteInitialData())),
               "initializeData()", "InitialData has to be written to preCICE before calling initializeData()");

  setHasDataBeenExchanged(false);

  if (hasToReceiveInitData() && isCouplingOngoing() )  {
    assertion(doesFirstStep());
    preciceDebug("Receiving data");
    startReceivePackage();
    receiveAndSetDt();
    receiveData(getCommunication());
    finishReceivePackage();
    setHasDataBeenExchanged(true);
  }

  if (hasToSendInitData() && isCouplingOngoing()) {
    assertion(not doesFirstStep());
    if (getExtrapolationOrder() > 0) {
      foreach (DataMap::value_type & pair, getSendData()) {
        if (pair.second->oldValues.cols() == 0)
          break;
        utils::DynVector& oldValues = pair.second->oldValues.column(0);
        oldValues = *pair.second->values;
        // For extrapolation, treat the initial value as old timestep value
        pair.second->oldValues.shiftSetFirst(*pair.second->values);
      }
    }
    // The second participant sends the initialized data to the first particpant
    // here, which receives the data on call of initialize().
    sendData(getCommunication());
    startReceivePackage();
    // This receive replaces the receive in initialize().
    receiveData(getCommunication());
    finishReceivePackage();
    setHasDataBeenExchanged(true);
  }

  //in order to check in advance if initializeData has been called (if necessary)
  setHasToSendInitData(false);
  setHasToReceiveInitData(false);
}

void SerialCouplingScheme:: advance()
{
  preciceTrace2("advance()", getTimesteps(), getTime());
  foreach (DataMap::value_type & pair, getReceiveData()) {
    utils::DynVector& values = *pair.second->values;
    preciceDebug("Begin advance, New Values: " << values);
  }
  checkCompletenessRequiredActions();

  preciceCheck(not hasToReceiveInitData() && not hasToSendInitData(), "advance()",
               "initializeData() needs to be called before advance if data has to be initialized!");
  
  setHasDataBeenExchanged(false);
  setIsCouplingTimestepComplete(false);

  if (_couplingMode == Explicit) {
    if (tarch::la::equals(getThisTimestepRemainder(), 0.0, _eps)) {
      setIsCouplingTimestepComplete(true);
      setTimesteps(getTimesteps() + 1);
      preciceDebug("Sending data...");
      startSendPackage();
      sendDt();
      sendData(getCommunication());
      finishSendPackage();

      if (isCouplingOngoing() || doesFirstStep()){
        preciceDebug("Receiving data...");
        startReceivePackage();
        receiveAndSetDt();
        receiveData(getCommunication());
        finishReceivePackage();
      }
      setHasDataBeenExchanged(true);
      setComputedTimestepPart(0.0);
    }
  }
  else if (_couplingMode == Implicit) {  
    bool convergence = true;
  
    if (tarch::la::equals(getThisTimestepRemainder(), 0.0, _eps)) {
      preciceDebug("Computed full length of iteration");
      if (doesFirstStep()) {
        startSendPackage();
        sendDt();
        sendData(getCommunication());
        finishSendPackage();
        startReceivePackage();
        getCommunication()->receive(convergence, 0); //TODO
        if (convergence) {
          timestepCompleted();
        }
        if (isCouplingOngoing()) {
          receiveData(getCommunication());
        }
        finishReceivePackage();
      }
      else {
        convergence = measureConvergence(); //TODO
        // Stop, when maximal iteration count (given in config) is reached
        if (maxIterationsReached()) {
          convergence = true;
        }
        if (convergence) {
          if (getPostProcessing().get() != NULL) {
            getPostProcessing()->iterationsConverged(getSendData());
          }
          newConvergenceMeasurements();
          timestepCompleted();
        }
        else if (getPostProcessing().get() != NULL) {
          getPostProcessing()->performPostProcessing(getSendData());
        }
        startSendPackage();
        getCommunication()->send(convergence, 0); //TODO
        if (isCouplingOngoing()) {
          if (convergence && (getExtrapolationOrder() > 0)){
            extrapolateData(getSendData()); // Also stores data
          }
          else { // Store data for conv. measurement, post-processing, or extrapolation
            foreach (DataMap::value_type& pair, getSendData()) {
              if (pair.second->oldValues.size() > 0){
                pair.second->oldValues.column(0) = *pair.second->values;
              }
            }
            foreach (DataMap::value_type& pair, getReceiveData()) {
              if (pair.second->oldValues.size() > 0){
                pair.second->oldValues.column(0) = *pair.second->values;
              }
            }
          }
          sendData(getCommunication());
          finishSendPackage();
          startReceivePackage();
          receiveAndSetDt();
          receiveData(getCommunication());
          finishReceivePackage();
        }
        else {
          finishSendPackage();
        }
      }
    
      if (not convergence) {
        preciceDebug("No convergence achieved");
        requireAction(constants::actionReadIterationCheckpoint());
      }
      else {
        preciceDebug("Convergence achieved");
        advanceTXTWriters();
      }
      updateTimeAndIterations(convergence);
      setHasDataBeenExchanged(true);
      setComputedTimestepPart(0.0);
    } //subcycling completed
    
  }
}



}}
