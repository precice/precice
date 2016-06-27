#include "MultiCouplingScheme.hpp"
#include "impl/PostProcessing.hpp"
#include "mesh/Mesh.hpp"
#include "com/Communication.hpp"
#include "utils/EigenHelperFunctions.hpp"
#include "m2n/M2N.hpp"

namespace precice {
namespace cplscheme {

logging::Logger MultiCouplingScheme::_log("precice::cplscheme::MultiCouplingScheme" );

MultiCouplingScheme::MultiCouplingScheme
(
  double                maxTime,
  int                   maxTimesteps,
  double                timestepLength,
  int                   validDigits,
  const std::string&    localParticipant,
  std::vector<m2n::M2N::SharedPointer> communications,
  constants::TimesteppingMethod dtMethod,
  int                   maxIterations)
  :
  BaseCouplingScheme(maxTime,maxTimesteps,timestepLength,validDigits,"neverFirstParticipant",
      localParticipant,localParticipant,m2n::M2N::SharedPointer(),maxIterations,dtMethod),
  _communications(communications),
  _allData (),
  _receiveDataVector(),
  _sendDataVector()
{
  for(size_t i = 0; i < _communications.size(); ++i) {
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
  preciceTrace("initialize()", startTime, startTimestep);
  assertion(not isInitialized());
  assertion(tarch::la::greaterEquals(startTime, 0.0), startTime);
  assertion(startTimestep >= 0, startTimestep);
  setTime(startTime);
  setTimesteps(startTimestep);


  mergeData(); // merge send and receive data for all pp calls
  setupConvergenceMeasures(); // needs _couplingData configured
  setupDataMatrices(_allData); // Reserve memory and initialize data with zero
  if (getPostProcessing().get() != nullptr) {
    preciceCheck(getPostProcessing()->getDataIDs().size()>=3 ,"initialize()",
                 "For parallel coupling, the number of coupling data vectors has to be at least 3, not: "
                 << getPostProcessing()->getDataIDs().size());
    getPostProcessing()->initialize(_allData); // Reserve memory, initialize
  }


  requireAction(constants::actionWriteIterationCheckpoint());
  initializeTXTWriters();


  for (DataMap& dataMap : _sendDataVector) {
    for (DataMap::value_type & pair : dataMap) {
      if (pair.second->initialize) {
        setHasToSendInitData(true);
        break;
      }
    }
  }
  for (DataMap& dataMap : _receiveDataVector) {
    for (DataMap::value_type & pair : dataMap) {
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
      for (DataMap& dataMap : _receiveDataVector) {
        for (DataMap::value_type & pair : dataMap){
          pair.second->oldValues.col(0) = *pair.second->values;
          // For extrapolation, treat the initial value as old timestep value
          utils::shiftSetFirst(pair.second->oldValues, *pair.second->values);
        }
      }
    }
  }
  if (hasToSendInitData()) {
    if (getExtrapolationOrder() > 0) {
      for (DataMap& dataMap : _sendDataVector) {
        for (DataMap::value_type & pair : dataMap) {
          pair.second->oldValues.col(0) = *pair.second->values;
          // For extrapolation, treat the initial value as old timestep value
          utils::shiftSetFirst(pair.second->oldValues, *pair.second->values);
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
  preciceTrace("advance()", getTimesteps(), getTime());
  checkCompletenessRequiredActions();

  preciceCheck(!hasToReceiveInitData() && !hasToSendInitData(), "advance()",
               "initializeData() needs to be called before advance if data has to be initialized!");

  setHasDataBeenExchanged(false);
  setIsCouplingTimestepComplete(false);
  bool convergence = false;
  if (tarch::la::equals(getThisTimestepRemainder(), 0.0, _eps)) {
    preciceDebug("Computed full length of iteration");

    receiveData();

    auto designSpecifications = getPostProcessing()->getDesignSpecification(_allData);
    convergence = measureConvergence(designSpecifications);

    // Stop, when maximal iteration count (given in config) is reached
    if (maxIterationsReached()) {
      convergence = true;
    }
    if (convergence) {
      if (getPostProcessing().get() != nullptr) {
        getPostProcessing()->iterationsConverged(_allData);
      }
      newConvergenceMeasurements();
      timestepCompleted();
    }
    else if (getPostProcessing().get() != nullptr) {
      getPostProcessing()->performPostProcessing(_allData);
    }

    for (m2n::M2N::SharedPointer m2n : _communications) {
      m2n->send(convergence);
      assertion(not _isCoarseModelOptimizationActive);
      m2n->send(_isCoarseModelOptimizationActive); //need to do this to match with ParallelCplScheme
    }

    if (convergence && (getExtrapolationOrder() > 0)){
      extrapolateData(_allData); // Also stores data
    }
    else { // Store data for conv. measurement, post-processing, or extrapolation
      for (DataMap::value_type& pair : _allData) {
        if (pair.second->oldValues.size() > 0){
          pair.second->oldValues.col(0) = *pair.second->values;
        }
      }
    }
    sendData();

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
  assertion(_allData.empty(), "This function should only be called once.");
  assertion(_sendDataVector.size()==_receiveDataVector.size());
  for(size_t i=0;i<_sendDataVector.size();i++){
    _allData.insert(_sendDataVector[i].begin(), _sendDataVector[i].end());
    _allData.insert(_receiveDataVector[i].begin(), _receiveDataVector[i].end());
  }
}

void MultiCouplingScheme:: addDataToSend
(
  mesh::PtrData data,
  mesh::PtrMesh mesh,
  bool          initialize,
  int           index)
{
  int id = data->getID();
  if(! utils::contained(id, _sendDataVector[index])) {
    PtrCouplingData ptrCplData (new CouplingData(& (data->values()), mesh, initialize, data->getDimensions()));
    DataMap::value_type pair = std::make_pair (id, ptrCplData);
    _sendDataVector[index].insert(pair);
  }
  else {
    preciceError("addDataToSend()", "Data \"" << data->getName()
     << "\" of mesh \"" << mesh->getName() << "\" cannot be "
     << "added twice for sending!");
  }
}

void MultiCouplingScheme:: addDataToReceive
(
  mesh::PtrData data,
  mesh::PtrMesh mesh,
  bool          initialize,
  int           index)
{
  int id = data->getID();
  if(! utils::contained(id, _receiveDataVector[index])) {
    PtrCouplingData ptrCplData (new CouplingData(& (data->values()), mesh, initialize, data->getDimensions()));
    DataMap::value_type pair = std::make_pair (id, ptrCplData);
    _receiveDataVector[index].insert(pair);
  }
  else {
    preciceError("addDataToReceive()", "Data \"" << data->getName()
     << "\" of mesh \"" << mesh->getName() << "\" cannot be "
     << "added twice for receiving!");
  }
}

void MultiCouplingScheme:: sendData()
{
  preciceTrace("sendData()");

  for(size_t i=0;i<_communications.size();i++){
    assertion(_communications[i].get() != nullptr);
    assertion(_communications[i]->isConnected());

    for (DataMap::value_type& pair : _sendDataVector[i]) {
      int size = pair.second->values->size();
      if (size > 0) {
        _communications[i]->send(pair.second->values->data(), size, pair.second->mesh->getID(), pair.second->dimension);
      }
    }
  }
}

void MultiCouplingScheme:: receiveData()
{
  preciceTrace("receiveData()");

  for(size_t i=0;i<_communications.size();i++){
    assertion(_communications[i].get() != nullptr);
    assertion(_communications[i]->isConnected());

    for (DataMap::value_type& pair : _receiveDataVector[i]) {
      int size = pair.second->values->size();
      if (size > 0) {
        _communications[i]->receive(pair.second->values->data(), size, pair.second->mesh->getID(), pair.second->dimension);
      }
    }
  }
}


void MultiCouplingScheme::setupConvergenceMeasures()
{
  preciceTrace("setupConvergenceMeasures()");
  assertion(not doesFirstStep());
  preciceCheck(not _convergenceMeasures.empty(), "setupConvergenceMeasures()",
         "At least one convergence measure has to be defined for "
         << "an implicit coupling scheme!");
  for (ConvergenceMeasure& convMeasure : _convergenceMeasures) {
    int dataID = convMeasure.dataID;
    convMeasure.data = getData(dataID);
    assertion(convMeasure.data != nullptr);
  }
}

CouplingData* MultiCouplingScheme:: getData
(
  int dataID)
{
  preciceTrace("getData()", dataID);
  DataMap::iterator iter = _allData.find(dataID);
  if (iter != _allData.end()) {
    return  &(*(iter->second));
  }
  return nullptr;
}

}}
