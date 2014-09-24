#include "MultiCouplingScheme.hpp"
#include "impl/PostProcessing.hpp"
#include "mesh/Mesh.hpp"
#include "com/Communication.hpp"
#include "com/SharedPointer.hpp"

namespace precice {
namespace cplscheme {

tarch::logging::Log MultiCouplingScheme::_log("precice::cplscheme::MultiCouplingScheme" );

MultiCouplingScheme::MultiCouplingScheme
(
  double                maxTime,
  int                   maxTimesteps,
  double                timestepLength,
  int                   validDigits,
  const std::string&    localParticipant,
  std::vector<com::PtrCommunication> communications,
  constants::TimesteppingMethod dtMethod,
  int                   maxIterations)
  :
  BaseCouplingScheme(maxTime,maxTimesteps,timestepLength,validDigits,"neverFirstParticipant",
      localParticipant,localParticipant,com::PtrCommunication(),maxIterations,dtMethod),
  _communications(communications),
  _allData (),
  _receiveDataVector(),
  _sendDataVector()
{
  for(size_t i=0;i<_communications.size();i++){
    DataMap receiveMap;
    DataMap sendMap;
    _receiveDataVector.push_back(receiveMap);
    _sendDataVector.push_back(sendMap);
  }
}

void MultiCouplingScheme::initialize
(
  double startTime,
  int    startTimestep )
{
  preciceTrace2("initialize()", startTime, startTimestep);
  assertion(not isInitialized());
  assertion1(tarch::la::greaterEquals(startTime, 0.0), startTime);
  assertion1(startTimestep >= 0, startTimestep);
  setTime(startTime);
  setTimesteps(startTimestep);

  setupConvergenceMeasures(); // needs _couplingData configured
  mergeData(); // merge send and receive data for all pp calls
  setupDataMatrices(_allData); // Reserve memory and initialize data with zero
  if (getPostProcessing().get() != NULL) {
    preciceCheck(getPostProcessing()->getDataIDs().size()>=3 ,"initialize()",
                 "For parallel coupling, the number of coupling data vectors has to be at least 3, not: "
                 << getPostProcessing()->getDataIDs().size());
    getPostProcessing()->initialize(_allData); // Reserve memory, initialize
  }


  requireAction(constants::actionWriteIterationCheckpoint());
  initializeTXTWriters();


  foreach (DataMap& dataMap, _sendDataVector){
    foreach (DataMap::value_type & pair, dataMap) {
      if (pair.second->initialize) {
        setHasToSendInitData(true);
        break;
      }
    }
  }
  foreach (DataMap& dataMap, _receiveDataVector){
    foreach (DataMap::value_type & pair, dataMap) {
      if (pair.second->initialize) {
        setHasToReceiveInitData(true);
        break;
      }
    }
  }

  if (hasToSendInitData()) {
    requireAction(constants::actionWriteInitialData());
  }

  setIsInitialized(true);
}

void MultiCouplingScheme::initializeData()
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


  if (hasToReceiveInitData()) {
    receiveData();
    setHasDataBeenExchanged(true);

    // second participant has to save values for extrapolation
    if (getExtrapolationOrder() > 0){
      foreach (DataMap& dataMap, _receiveDataVector){
        foreach (DataMap::value_type & pair, dataMap){
          utils::DynVector& oldValues = pair.second->oldValues.column(0);
          oldValues = *pair.second->values;
          // For extrapolation, treat the initial value as old timestep value
          pair.second->oldValues.shiftSetFirst(*pair.second->values);
        }
      }
    }
  }
  if (hasToSendInitData()) {
    if (getExtrapolationOrder() > 0) {
      foreach (DataMap& dataMap, _sendDataVector){
        foreach (DataMap::value_type & pair, dataMap) {
          utils::DynVector& oldValues = pair.second->oldValues.column(0);
          oldValues = *pair.second->values;
          // For extrapolation, treat the initial value as old timestep value
          pair.second->oldValues.shiftSetFirst(*pair.second->values);
        }
      }
    }
    sendData();
  }


  // in order to check in advance if initializeData has been called (if necessary)
  setHasToSendInitData(false);
  setHasToReceiveInitData(false);
}

void MultiCouplingScheme::advance()
{
  preciceTrace2("advance()", getTimesteps(), getTime());
  checkCompletenessRequiredActions();

  preciceCheck(!hasToReceiveInitData() && !hasToSendInitData(), "advance()",
               "initializeData() needs to be called before advance if data has to be initialized!");

  setHasDataBeenExchanged(false);
  setIsCouplingTimestepComplete(false);
  bool convergence = false;
  if (tarch::la::equals(getThisTimestepRemainder(), 0.0, _eps)) {
    preciceDebug("Computed full length of iteration");

    receiveData();

    convergence = measureConvergence();

    // Stop, when maximal iteration count (given in config) is reached
    if (maxIterationsReached()) {
      convergence = true;
    }
    if (convergence) {
      if (getPostProcessing().get() != NULL) {
        getPostProcessing()->iterationsConverged(_allData);
      }
      newConvergenceMeasurements();
      timestepCompleted();
    }
    else if (getPostProcessing().get() != NULL) {
      getPostProcessing()->performPostProcessing(_allData);
    }

    foreach(com::PtrCommunication com, _communications){
      com->send(convergence, 0);
    }

    if (isCouplingOngoing()) {
      if (convergence && (getExtrapolationOrder() > 0)){
        extrapolateData(_allData); // Also stores data
      }
      else { // Store data for conv. measurement, post-processing, or extrapolation
        foreach (DataMap::value_type& pair, _allData) {
          if (pair.second->oldValues.size() > 0){
            pair.second->oldValues.column(0) = *pair.second->values;
          }
        }
      }
      sendData();
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
  } // subcycling complete
}




void MultiCouplingScheme::mergeData()
{
  preciceTrace("mergeData()");
  assertion1(_allData.empty(), "This function should only be called once.");
  assertion(_sendDataVector.size()==_receiveDataVector.size());
  for(size_t i=0;i<_sendDataVector.size();i++){
    _allData.insert(_sendDataVector[i].begin(), _sendDataVector[i].end());
    _allData.insert(_receiveDataVector[i].begin(), _receiveDataVector[i].end());
  }
}

void MultiCouplingScheme:: addDataToSend
(
  mesh::PtrData data,
  bool          initialize,
  int           index)
{
  int id = data->getID();
  if(! utils::contained(id, _sendDataVector[index])) {
    PtrCouplingData ptrCplData (new CouplingData(& (data->values()), initialize));
    DataMap::value_type pair = std::make_pair (id, ptrCplData);
    _sendDataVector[index].insert(pair);
  }
  else {
    preciceError("addDataToSend()", "Data \"" << data->getName()
     << "\" of mesh \"" << (data->mesh())->getName() << "\" cannot be "
     << "added twice for sending!");
  }
}

void MultiCouplingScheme:: addDataToReceive
(
  mesh::PtrData data,
  bool          initialize,
  int           index)
{
  int id = data->getID();
  if(! utils::contained(id, _receiveDataVector[index])) {
    PtrCouplingData ptrCplData (new CouplingData(& (data->values()), initialize));
    DataMap::value_type pair = std::make_pair (id, ptrCplData);
    _receiveDataVector[index].insert(pair);
  }
  else {
    preciceError("addDataToReceive()", "Data \"" << data->getName()
     << "\" of mesh \"" << (data->mesh())->getName() << "\" cannot be "
     << "added twice for receiving!");
  }
}

void MultiCouplingScheme:: sendData()
{
  preciceTrace("sendData()");

  for(size_t i=0;i<_communications.size();i++){
    assertion(_communications[i].get() != NULL);
    assertion(_communications[i]->isConnected());

    foreach (DataMap::value_type& pair, _sendDataVector[i]){
      int size = pair.second->values->size();
      if (size > 0) {
        _communications[i]->send(tarch::la::raw(*pair.second->values), size, 0);
      }
    }
  }
}

void MultiCouplingScheme:: receiveData()
{
  preciceTrace("receiveData()");

  for(size_t i=0;i<_communications.size();i++){
    assertion(_communications[i].get() != NULL);
    assertion(_communications[i]->isConnected());

    foreach (DataMap::value_type& pair, _receiveDataVector[i]){
      int size = pair.second->values->size();
      if (size > 0) {
        _communications[i]->receive(tarch::la::raw(*pair.second->values), size, 0);
      }
    }
  }
}

}}
