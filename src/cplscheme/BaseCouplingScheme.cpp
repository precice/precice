// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "BaseCouplingScheme.hpp"
#include "mesh/Mesh.hpp"
#include "com/Communication.hpp"
#include "m2n/M2N.hpp"
#include "utils/Globals.hpp"
#include "utils/MasterSlave.hpp"
#include "impl/PostProcessing.hpp"
#include "impl/ConvergenceMeasure.hpp"
#include "io/TXTWriter.hpp"
#include "io/TXTReader.hpp"
#include <limits>
#include <sstream>

namespace precice {
namespace cplscheme {

tarch::logging::Log BaseCouplingScheme::
_log("precice::cplscheme::BaseCouplingScheme");

BaseCouplingScheme:: BaseCouplingScheme
(
  double maxTime,
  int    maxTimesteps,
  double timestepLength,
  int    validDigits)
  :
  _couplingMode(Undefined),
  _eps(std::pow(10.0, -1 * validDigits)),
  _participantSetsDt(false),
  _participantReceivesDt(false),
  _maxTime(maxTime),
  _maxTimesteps(maxTimesteps),
  _iterations(-1),
  _maxIterations(-1),
  _totalIterations(-1),
  _timesteps(0),
  _timestepLength(timestepLength),
  _time(0.0),
  _computedTimestepPart(0.0),
  _extrapolationOrder(0),
  _validDigits(validDigits),
  _firstResiduumNorm(0),
  _doesFirstStep(false),
  _checkpointTimestepInterval(-1),
  _isCouplingTimestepComplete(false),
  _hasToSendInitData(false),
  _hasToReceiveInitData(false),
  _hasDataBeenExchanged(false),
  _isInitialized(false),
  _actions(),
  _sendData(),
  _receiveData (),
  _iterationsWriter("iterations-unknown.txt"),
  _convergenceWriter("convergence-unknown.txt")
{
  preciceCheck (
    not ((maxTime != UNDEFINED_TIME) && (maxTime < 0.0)),
    "BaseCouplingScheme()", "Maximum time has to be larger than zero!");
  preciceCheck (
    not ((maxTimesteps != UNDEFINED_TIMESTEPS) && (maxTimesteps < 0)),
    "BaseCouplingScheme()", "Maximum timestep number has to be larger than zero!");
  preciceCheck (
    not ((timestepLength != UNDEFINED_TIMESTEP_LENGTH) && (timestepLength < 0.0)),
    "BaseCouplingScheme()", "Timestep length has to be larger than zero!");
  preciceCheck((_validDigits >= 1) && (_validDigits < 17),
           "BaseCouplingScheme()", "Valid digits of timestep length has to be "
           << "between 1 and 16!");
}

BaseCouplingScheme::BaseCouplingScheme
(
  double                maxTime,
  int                   maxTimesteps,
  double                timestepLength,
  int                   validDigits,
  const std::string&    firstParticipant,
  const std::string&    secondParticipant,
  const std::string&    localParticipant,
  m2n::M2N::SharedPointer           m2n,
  int                   maxIterations,
  constants::TimesteppingMethod dtMethod )
  :
  _firstParticipant(firstParticipant),
  _secondParticipant(secondParticipant),
  _convergenceMeasures(),
  _eps(std::pow(10.0, -1 * validDigits)),
  _m2n(m2n),
  _participantSetsDt(false),
  _participantReceivesDt(false),
  _maxTime(maxTime),
  _maxTimesteps(maxTimesteps),
  _iterations(1),
  _maxIterations(maxIterations),
  _totalIterations(1),
  _timesteps(1),
  _timestepLength(timestepLength),
  _time(0.0),
  _computedTimestepPart(0.0),
  _extrapolationOrder(0),
  _firstResiduumNorm(0),
  _validDigits(validDigits),
  _doesFirstStep(false),
  _checkpointTimestepInterval(-1),
  _isCouplingTimestepComplete(false),
  _postProcessing(),
  _hasToSendInitData(false),
  _hasToReceiveInitData(false),
  _hasDataBeenExchanged(false),
  _isInitialized(false),
  _actions(),
  _sendData(),
  _receiveData(),
  _iterationsWriter("iterations-" + localParticipant + ".txt"),
  _convergenceWriter("convergence-" + localParticipant + ".txt")
{
  preciceCheck(
    not ((maxTime != UNDEFINED_TIME) && (maxTime < 0.0)),
    "BaseCouplingScheme()", "Maximum time has to be larger than zero!");
  preciceCheck(
    not ((maxTimesteps != UNDEFINED_TIMESTEPS) && (maxTimesteps < 0)),
    "BaseCouplingScheme()", "Maximum timestep number has to be larger than zero!");
  preciceCheck(
    not ((timestepLength != UNDEFINED_TIMESTEP_LENGTH) && (timestepLength < 0.0)),
    "BaseCouplingScheme()", "Timestep length has to be larger than zero!");
  preciceCheck((_validDigits >= 1) && (_validDigits < 17),
           "BaseCouplingScheme()", "Valid digits of timestep length has to be "
               << "between 1 and 16!");
  preciceCheck(_firstParticipant != _secondParticipant,
           "ImplicitCouplingScheme()", "First participant and "
           << "second participant must have different names! Called from BaseCoupling.");
  if (dtMethod == constants::FIXED_DT){
    preciceCheck(hasTimestepLength(), "ImplicitCouplingScheme()",
         "Timestep length value has to be given "
         << "when the fixed timestep length method is chosen for an implicit "
         << "coupling scheme!");
  }
  if (localParticipant == _firstParticipant){
    _doesFirstStep = true;
    if (dtMethod == constants::FIRST_PARTICIPANT_SETS_DT){
      _participantSetsDt = true;
      setTimestepLength(UNDEFINED_TIMESTEP_LENGTH);
    }
  }
  else if (localParticipant == _secondParticipant){
    if (dtMethod == constants::FIRST_PARTICIPANT_SETS_DT){
      _participantReceivesDt = true;
    }
  }
  else {
    preciceError("initialize()", "Name of local participant \""
         << localParticipant << "\" does not match any "
         << "participant specified for the coupling scheme!");
  }
  preciceCheck((maxIterations > 0) || (maxIterations == -1),
           "ImplicitCouplingState()",
           "Maximal iteration limit has to be larger than zero!");
}


void BaseCouplingScheme:: receiveAndSetDt()
{
  preciceTrace("receiveAndSetDt()");
  if (participantReceivesDt()){
    double dt = UNDEFINED_TIMESTEP_LENGTH;
    getM2N()->receive(dt);
    preciceDebug("Received timestep length of " << dt);
    assertion(not tarch::la::equals(dt, UNDEFINED_TIMESTEP_LENGTH));
    setTimestepLength(dt);
  }
}

void BaseCouplingScheme:: sendDt(){
  preciceTrace("sendDt()");
  if (participantSetsDt()){
    preciceDebug("sending timestep length of " << getComputedTimestepPart());
    getM2N()->send(getComputedTimestepPart());
  }
}


void BaseCouplingScheme:: addDataToSend
(
  mesh::PtrData data,
  mesh::PtrMesh mesh,
  bool          initialize)
{
  preciceTrace("addDataToSend()");
  int id = data->getID();
  if(! utils::contained(id, _sendData)) {
    PtrCouplingData ptrCplData (new CouplingData(& (data->values()), mesh, initialize, data->getDimensions()));
    DataMap::value_type pair = std::make_pair (id, ptrCplData);
    _sendData.insert(pair);
  }
  else {
    preciceError("addDataToSend()", "Data \"" << data->getName()
         << "\" cannot be added twice for sending!");
  }
}

void BaseCouplingScheme:: addDataToReceive
(
  mesh::PtrData data,
  mesh::PtrMesh mesh,
  bool          initialize)
{
  preciceTrace("addDataToReceive()");
  int id = data->getID();
  if(! utils::contained(id, _receiveData)) {
    PtrCouplingData ptrCplData (new CouplingData(& (data->values()), mesh, initialize, data->getDimensions()));
    DataMap::value_type pair = std::make_pair (id, ptrCplData);
    _receiveData.insert(pair);
  }
  else {
    preciceError("addDataToReceive()", "Data \"" << data->getName()
         << "\" cannot be added twice for receiving!");
  }
}


void BaseCouplingScheme:: sendState
(
  com::Communication::SharedPointer communication,
  int                   rankReceiver)
{
  preciceTrace1("sendState()", rankReceiver);
  communication->startSendPackage(rankReceiver );
  assertion(communication.get() != NULL);
  assertion(communication->isConnected());
  communication->send(_maxTime, rankReceiver);
  communication->send(_maxTimesteps, rankReceiver);
  communication->send(_timestepLength, rankReceiver);
  communication->send(_time, rankReceiver);
  communication->send(_timesteps, rankReceiver);
  communication->send(_checkpointTimestepInterval, rankReceiver);
  communication->send(_computedTimestepPart, rankReceiver);
  //communication->send(_maxLengthNextTimestep, rankReceiver);
  communication->send(_isInitialized, rankReceiver);
  communication->send(_isCouplingTimestepComplete, rankReceiver);
  communication->send(_hasDataBeenExchanged, rankReceiver);
  communication->send((int)_actions.size(), rankReceiver);
  for (const std::string& action : _actions) {
    communication->send(action, rankReceiver);
  }
  communication->send(_maxIterations, rankReceiver );
  communication->send(_iterations, rankReceiver );
  communication->send(_totalIterations, rankReceiver );
  communication->finishSendPackage();

}

void BaseCouplingScheme:: receiveState
(
  com::Communication::SharedPointer communication,
  int                   rankSender)
{
  preciceTrace1("receiveState()", rankSender);
  communication->startReceivePackage(rankSender);
  assertion(communication.get() != NULL);
  assertion(communication->isConnected());
  communication->receive(_maxTime, rankSender);
  communication->receive(_maxTimesteps, rankSender);
  communication->receive(_timestepLength, rankSender);
  communication->receive(_time, rankSender);
  communication->receive(_timesteps, rankSender);
  communication->receive(_checkpointTimestepInterval, rankSender);
  communication->receive(_computedTimestepPart, rankSender);
  //communication->receive(_maxLengthNextTimestep, rankSender);
  communication->receive(_isInitialized, rankSender);
  communication->receive(_isCouplingTimestepComplete, rankSender);
  communication->receive(_hasDataBeenExchanged, rankSender);
  int actionsSize = 0;
  communication->receive(actionsSize, rankSender);
  _actions.clear();
  for (int i=0; i < actionsSize; i++) {
    std::string action;
    communication->receive(action, rankSender);
    _actions.insert(action);
  }
  communication->receive(_maxIterations, rankSender);
  int subIteration = -1;
  communication->receive(subIteration, rankSender);
  _iterations = subIteration;
  communication->receive(_totalIterations, rankSender);
  communication->finishReceivePackage();

}

std::vector<int> BaseCouplingScheme:: sendData
(
  m2n::M2N::SharedPointer m2n)
{
  preciceTrace("sendData()");

  std::vector<int> sentDataIDs;
  assertion(m2n.get() != NULL);
  assertion(m2n->isConnected());
  for (DataMap::value_type& pair : _sendData){
    int size = pair.second->values->size();
    m2n->send(tarch::la::raw(*(pair.second->values)), size,
              pair.second->mesh->getID(), pair.second->dimension);
    sentDataIDs.push_back(pair.first);
  }
  preciceDebug("Number of sent data sets = " << sentDataIDs.size());
  return sentDataIDs;
}

std::vector<int> BaseCouplingScheme:: receiveData
(
  m2n::M2N::SharedPointer m2n)
{
  preciceTrace("receiveData()");
  std::vector<int> receivedDataIDs;
  assertion(m2n.get() != NULL);
  assertion(m2n->isConnected());

  for (DataMap::value_type & pair : _receiveData) {
    int size = pair.second->values->size();
    m2n->receive(tarch::la::raw(*(pair.second->values)), size,
                 pair.second->mesh->getID(), pair.second->dimension);
    receivedDataIDs.push_back(pair.first);
  }
  preciceDebug("Number of received data sets = " << receivedDataIDs.size());

  return receivedDataIDs;
}


int BaseCouplingScheme:: getVertexOffset(
    std::map<int,int>& vertexDistribution,
    int rank,
    int dim)
  {
    int sum=0;
    for(int i=0;i<rank;i++){
      sum += vertexDistribution[i];
    }
    return sum*dim;
  }

CouplingData* BaseCouplingScheme:: getSendData
(
  int dataID)
{
  preciceTrace1("getSendData()", dataID);
  DataMap::iterator iter = _sendData.find(dataID);
  if (iter != _sendData.end()) {
    return  &(*(iter->second));
  }
  return NULL;
}

CouplingData* BaseCouplingScheme:: getReceiveData
(
  int dataID)
{
  preciceTrace1("getReceiveData()", dataID);
  DataMap::iterator iter = _receiveData.find(dataID);
  if (iter != _receiveData.end()) {
    return  &(*(iter->second));
  }
  return NULL;
}

void BaseCouplingScheme::finalize()
{
  preciceTrace("finalize()");
  checkCompletenessRequiredActions();
  preciceCheck(isInitialized(), "finalize()",
           "Called finalize() before initialize()!");
}

void BaseCouplingScheme:: setExtrapolationOrder
(
  int order)
{
  preciceCheck((order == 0) || (order == 1) || (order == 2),
               "setExtrapolationOrder()", "Extrapolation order has to be "
               << " 0, 1, or 2!");
  _extrapolationOrder = order;
}


void BaseCouplingScheme::extrapolateData(DataMap& data)
{
  preciceTrace1("extrapolateData()", _timesteps);
  if ((_extrapolationOrder == 1) || getTimesteps() == 2) { //timesteps is increased before extrapolate is called
    preciceInfo("extrapolateData()", "Performing first order extrapolation" );
    for (DataMap::value_type & pair : data) {
      preciceDebug("Extrapolate data: " << pair.first);
      assertion(pair.second->oldValues.cols() > 1 );
      utils::DynVector & values = *pair.second->values;
      pair.second->oldValues.column(0) = values;     // = x^t
      values *= 2.0;                                 // = 2*x^t
      values -= pair.second->oldValues.column(1);    // = 2*x^t - x^(t-1)
      pair.second->oldValues.shiftSetFirst(values ); // shift old values to the right
    }
  }
  else if (_extrapolationOrder == 2 ) {
    preciceInfo("extrapolateData()", "Performing second order extrapolation" );
    for (DataMap::value_type & pair : data ) {
      assertion(pair.second->oldValues.cols() > 2 );
      utils::DynVector & values = *pair.second->values;
      utils::DynVector & valuesOld1 = pair.second->oldValues.column(1);
      utils::DynVector & valuesOld2 = pair.second->oldValues.column(2);

      pair.second->oldValues.column(0) = values;        // = x^t
      values *= 2.5;                                    // = 2.5 x^t
      values -= valuesOld1 * 2.0; // = 2.5x^t - 2x^(t-1)
      values += valuesOld2 * 0.5; // = 2.5x^t - 2x^(t-1) + 0.5x^(t-2)
      pair.second->oldValues.shiftSetFirst(values);
    }
  }
  else {
    preciceError("extrapolateData()", "Called extrapolation with order != 1,2!" );
  }
}

bool BaseCouplingScheme:: hasTimestepLength() const
{
  return not tarch::la::equals(_timestepLength, UNDEFINED_TIMESTEP_LENGTH);
}

double BaseCouplingScheme:: getTimestepLength() const
{
  assertion(not tarch::la::equals(_timestepLength, UNDEFINED_TIMESTEP_LENGTH));
  return _timestepLength;
}

void BaseCouplingScheme:: addComputedTime
(
  double timeToAdd )
{
  preciceTrace2("addComputedTime()", timeToAdd, _time);
  preciceCheck(isCouplingOngoing(), "addComputedTime()",
           "Invalid call of addComputedTime() after simulation end!");

  _computedTimestepPart += timeToAdd;
  _time += timeToAdd;

  // Check validness
  bool valid = tarch::la::greaterEquals(getThisTimestepRemainder(), 0.0, _eps);
  preciceCheck(valid, "addComputedTime()", "The computed timestep length of "
           << timeToAdd << " exceeds the maximum timestep limit of "
           << _timestepLength - _computedTimestepPart + timeToAdd
           << " for this time step!");
}

bool BaseCouplingScheme:: willDataBeExchanged
(
  double lastSolverTimestepLength) const
{
  preciceTrace1("willDataBeExchanged()", lastSolverTimestepLength);
  double remainder = getThisTimestepRemainder() - lastSolverTimestepLength;
  return not tarch::la::greater(remainder, 0.0, _eps);
}

bool BaseCouplingScheme:: hasDataBeenExchanged() const
{
  return _hasDataBeenExchanged;
}

void BaseCouplingScheme:: setHasDataBeenExchanged
(
  bool hasDataBeenExchanged)
{
  _hasDataBeenExchanged = hasDataBeenExchanged;
}

double BaseCouplingScheme:: getTime() const
{
  return _time;
}

int BaseCouplingScheme:: getTimesteps() const
{
  return _timesteps;
}


std::vector<std::string> BaseCouplingScheme::getCouplingPartners() const
{
  std::vector<std::string> partnerNames;
  // Add non-local participant
  if (doesFirstStep()) {
    partnerNames.push_back(_secondParticipant);
  }
  else {
    partnerNames.push_back(_firstParticipant);
  }
  return partnerNames;
}


double BaseCouplingScheme:: getThisTimestepRemainder() const
{
  preciceTrace("getTimestepRemainder()");
  double remainder = 0.0;
  if (not tarch::la::equals(_timestepLength, UNDEFINED_TIMESTEP_LENGTH)){
    remainder = _timestepLength - _computedTimestepPart;
  }
  preciceDebug("return " << remainder);
  return remainder;
}

double BaseCouplingScheme:: getNextTimestepMaxLength() const
{
  if (tarch::la::equals(_timestepLength, UNDEFINED_TIMESTEP_LENGTH)){
    if (tarch::la::equals(_maxTime, UNDEFINED_TIME)){
      return std::numeric_limits<double>::max();
    }
    else {
      return _maxTime - _time;
    }
  }
  return _timestepLength - _computedTimestepPart;
}

bool BaseCouplingScheme:: isCouplingOngoing() const
{
  using namespace tarch::la;
  bool timeLeft = greater(_maxTime, _time, _eps) || equals(_maxTime, UNDEFINED_TIME);
  bool timestepsLeft = (_maxTimesteps >= _timesteps)
    || (_maxTimesteps == UNDEFINED_TIMESTEPS);
  return timeLeft && timestepsLeft;
}

bool BaseCouplingScheme:: isCouplingTimestepComplete() const
{
  return _isCouplingTimestepComplete;
}

bool BaseCouplingScheme:: isActionRequired
(
  const std::string& actionName) const
{
  return _actions.count(actionName) > 0;
}

void BaseCouplingScheme:: performedAction
(
  const std::string& actionName)
{
  _actions.erase(actionName);
}

int BaseCouplingScheme:: getCheckpointTimestepInterval() const
{
  return _checkpointTimestepInterval;
}

void BaseCouplingScheme:: requireAction
(
  const std::string& actionName)
{
  _actions.insert(actionName);
}

std::string BaseCouplingScheme::printCouplingState() const
{
  std::ostringstream os;
  os << "it " << _iterations; //_iterations;
  if (getMaxIterations() != -1 ) {
    os << " of " << getMaxIterations();
  }
  os << " | " << printBasicState(_timesteps, _time) << " | " << printActionsState();
  return os.str();
}

std::string BaseCouplingScheme:: printBasicState() const
{
  std::ostringstream os;
  os << printBasicState(_timesteps, _time);
  return os.str ();
}

std::string BaseCouplingScheme:: printBasicState
(
  int    timesteps,
  double time ) const
{
  std::ostringstream os;
  os << "dt# " << timesteps;
  if(_maxTimesteps != UNDEFINED_TIMESTEPS){
    os << " of " << _maxTimesteps;
  }
  os << " | t " << time;
  if(_maxTime != UNDEFINED_TIME){
    os << " of " << _maxTime;
  }
  if(_timestepLength != UNDEFINED_TIMESTEP_LENGTH){
    os << " | dt " << _timestepLength;
  }
  if((_timestepLength != UNDEFINED_TIMESTEP_LENGTH)
     || (_maxTime != UNDEFINED_TIME))
  {
    os << " | max dt " << getNextTimestepMaxLength();
  }
  os << " | ongoing ";
  isCouplingOngoing() ? os << "yes" : os << "no";
  os << " | dt complete ";
  _isCouplingTimestepComplete ? os << "yes" : os << "no";
  return os.str ();
}

std::string BaseCouplingScheme:: printActionsState () const
{
  std::ostringstream os;
  for (const std::string & actionName : _actions) {
    os << actionName << " | ";
  }
  return os.str ();
}

void BaseCouplingScheme:: checkCompletenessRequiredActions ()
{
  preciceTrace("checkCompletenessRequiredActions()");
  if(not _actions.empty()){
    std::ostringstream stream;
    for (const std::string & action : _actions) {
      if (not stream.str().empty()){
    stream << ", ";
      }
      stream << action;
    }
    preciceError("checkCompletenessRequiredActions()",
         "Unfulfilled required actions: " << stream.str() << "!");
  }
}

int BaseCouplingScheme:: getValidDigits () const
{
  return _validDigits;
}

void BaseCouplingScheme::setupDataMatrices(DataMap& data)
{
  preciceTrace("setupDataMatrices()");
  preciceDebug("Data size: " << data.size());
  // Reserve storage for convergence measurement of send and receive data values
  for (ConvergenceMeasure& convMeasure : _convergenceMeasures) {
    assertion(convMeasure.data != NULL);
    if (convMeasure.data->oldValues.cols() < 1){
      convMeasure.data->oldValues.append(CouplingData::DataMatrix(
                       convMeasure.data->values->size(), 1, 0.0));
    }
  }
  // Reserve storage for extrapolation of data values
  if (_extrapolationOrder > 0){
    for (DataMap::value_type& pair : data) {
      int cols = pair.second->oldValues.cols();
      preciceDebug("Add cols: " << pair.first << ", cols: " << cols);
      assertion1(cols <= 1, cols);
      pair.second->oldValues.append(CouplingData::DataMatrix(
                      pair.second->values->size(), _extrapolationOrder + 1 - cols, 0.0));
    }
  }
}

void BaseCouplingScheme::setIterationPostProcessing
(
  impl::PtrPostProcessing postProcessing )
{
  assertion(postProcessing.get() != NULL);
  _postProcessing = postProcessing;
}


void BaseCouplingScheme::setupConvergenceMeasures()
{
  preciceTrace("setupConvergenceMeasures()");
  assertion(not doesFirstStep());
  preciceCheck(not _convergenceMeasures.empty(), "setupConvergenceMeasures()",
           "At least one convergence measure has to be defined for "
           << "an implicit coupling scheme!");
  for (ConvergenceMeasure& convMeasure : _convergenceMeasures) {
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

void BaseCouplingScheme::newConvergenceMeasurements()
{
  preciceTrace("newConvergenceMeasurements()");
  for (ConvergenceMeasure& convMeasure : _convergenceMeasures) {
    assertion(convMeasure.measure.get() != NULL);
    convMeasure.measure->newMeasurementSeries();
  }
}


void BaseCouplingScheme::addConvergenceMeasure
(
  int                         dataID,
  bool                        suffices,
  impl::PtrConvergenceMeasure measure )
{
  ConvergenceMeasure convMeasure;
  convMeasure.dataID = dataID;
  convMeasure.data = NULL;
  convMeasure.suffices = suffices;
  convMeasure.measure = measure;
  _convergenceMeasures.push_back(convMeasure);
  _firstResiduumNorm.push_back(0);
}


bool BaseCouplingScheme:: measureConvergence()
{
  preciceTrace("measureConvergence()");
  bool allConverged = true;
  bool oneSuffices = false;
  assertion(_convergenceMeasures.size() > 0);
  _convergenceWriter.writeData("Timestep", _timesteps);
  _convergenceWriter.writeData("Iteration", _iterations);
  for(int i = 0; i < _convergenceMeasures.size(); i++) {
    ConvergenceMeasure& convMeasure = _convergenceMeasures[i];
//  for (ConvergenceMeasure& convMeasure : _convergenceMeasures) {
    assertion(convMeasure.data != NULL);
    assertion(convMeasure.measure.get() != NULL);
    utils::DynVector& oldValues = convMeasure.data->oldValues.column(0);
    convMeasure.measure->measure(oldValues, *convMeasure.data->values);
    std::stringstream sstm;
    sstm << "resNorm(" <<i<< ")";
    _convergenceWriter.writeData(sstm.str(), convMeasure.measure->getNormResidual());
    if(_iterations == 1) {
      _firstResiduumNorm[i] = convMeasure.measure->getNormResidual(); 
    }
    if (not convMeasure.measure->isConvergence()) {
      //preciceDebug("Local convergence = false");
      allConverged = false;
    }
    else if (convMeasure.suffices == true) {
      oneSuffices = true;
    }
    preciceInfo("measureConvergence()", convMeasure.measure->printState());
  }
  if (allConverged) {
    preciceInfo("measureConvergence()", "All converged");
  }
  else if (oneSuffices) {
    preciceInfo("measureConvergence()", "Sufficient measure converged");
  }
  return allConverged || oneSuffices;
}



void BaseCouplingScheme::initializeTXTWriters()
{
  if(not utils::MasterSlave::_slaveMode){
    _iterationsWriter.addData("Timesteps", io::TXTTableWriter::INT );
    _iterationsWriter.addData("Total Iterations", io::TXTTableWriter::INT );
    _iterationsWriter.addData("Iterations", io::TXTTableWriter::INT );
    _iterationsWriter.addData("Convergence", io::TXTTableWriter::INT );
    int i = 0;
    for (ConvergenceMeasure& convMeasure : _convergenceMeasures) 
    {
       std::stringstream sstm; 
       sstm << "avgConvRate(" <<i<<")";
       _iterationsWriter.addData(sstm.str(), io::TXTTableWriter::DOUBLE);
       i++;
    }

    _convergenceWriter.addData("Timestep", io::TXTTableWriter::INT );
    _convergenceWriter.addData("Iteration", io::TXTTableWriter::INT );
     i = 0;
    for (ConvergenceMeasure& convMeasure : _convergenceMeasures) 
    {
       std::stringstream sstm;
       sstm << "resNorm(" <<i<< ")";
       _convergenceWriter.addData(sstm.str(), io::TXTTableWriter::DOUBLE);                                                                                                                                                         
       i++;
    }
  }
}

void BaseCouplingScheme::advanceTXTWriters()
{
  if(not utils::MasterSlave::_slaveMode){
    _iterationsWriter.writeData("Timesteps", _timesteps-1);
    _iterationsWriter.writeData("Total Iterations", _totalIterations);
    _iterationsWriter.writeData("Iterations", _iterations);
    int converged = _iterations < _maxIterations ? 1 : 0;
    _iterationsWriter.writeData("Convergence", converged);
    for (int i = 0; i<_convergenceMeasures.size();i++) {
      double avgConvRate = _convergenceMeasures[i].measure->getNormResidual()/_firstResiduumNorm[i];
      std::stringstream sstm;
      sstm << "avgConvRate(" <<i<< ")";
      _iterationsWriter.writeData(sstm.str(), std::pow(avgConvRate, 1./(double)_iterations));
    }
  }
}


void BaseCouplingScheme:: exportState(const std::string& filenamePrefix ) const
{
  if (not doesFirstStep()) {
    io::TXTWriter writer(filenamePrefix + "_cplscheme.txt");
    for (const BaseCouplingScheme::DataMap::value_type& dataMap : getSendData()) {
      writer.write(dataMap.second->oldValues);
    }
    for (const BaseCouplingScheme::DataMap::value_type& dataMap : getReceiveData()) {
      writer.write(dataMap.second->oldValues);
    }
    if (_postProcessing.get() != NULL) {
      _postProcessing->exportState(writer);
    }
  }
}

void BaseCouplingScheme:: importState(const std::string& filenamePrefix)
{
  if (not doesFirstStep()) {
    io::TXTReader reader(filenamePrefix + "_cplscheme.txt");
    for (BaseCouplingScheme::DataMap::value_type& dataMap : getSendData()) {
      reader.read(dataMap.second->oldValues);
    }
    for (BaseCouplingScheme::DataMap::value_type& dataMap : getReceiveData()) {
      reader.read(dataMap.second->oldValues);
    }
    if (_postProcessing.get() != NULL){
      _postProcessing->importState(reader);
    }
  }
}


void BaseCouplingScheme:: updateTimeAndIterations(bool convergence){
  _totalIterations++;
  if(not convergence){
    _iterations++;
    // The computed timestep part equals the timestep length, since the
    // timestep remainder is zero. Subtract the timestep length do another
    // coupling iteration.
    assertion(tarch::la::greater(getComputedTimestepPart(), 0.0));
    _time = _time - _computedTimestepPart;
  } else{
    _iterations = 1;
  }
}

void BaseCouplingScheme:: timestepCompleted()
{
  preciceTrace2("timestepCompleted()", getTimesteps(), getTime());
  preciceInfo("timestepCompleted()", "Timestep completed");
  setIsCouplingTimestepComplete(true);
  setTimesteps(getTimesteps() + 1 );
  //setTime(getTimesteps() * getTimestepLength() ); // Removes numerical errors
  if (isCouplingOngoing()) {
    preciceDebug("Setting require create checkpoint");
    requireAction(constants::actionWriteIterationCheckpoint());
  }
}

bool BaseCouplingScheme:: maxIterationsReached(){
  return _iterations == _maxIterations;
}




}} // namespace precice, cplscheme
