// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "BaseCouplingScheme.hpp"
#include "CompositionalCouplingScheme.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Data.hpp"
#include "com/Communication.hpp"
#include "utils/Globals.hpp"
#include "utils/Parallel.hpp"
#include <set>
#include <limits>

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
  _maxTime(maxTime),
  _maxTimesteps(maxTimesteps),
  _timestepLength(timestepLength),
  _validDigits(validDigits),
  _time(0.0),
  _computedTimestepPart(0.0),
  _timesteps(0),
  _checkpointTimestepInterval(-1),
  _isCouplingOngoing(true),
  _isCouplingTimestepComplete(false),
  _hasDataBeenExchanged(false),
  _isInitialized(false),
  _actions (),
  _sendData (),
  _receiveData ()
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

void BaseCouplingScheme:: addDataToSend
(
  mesh::PtrData data,
  bool          initialize)
{
  int id = data->getID();
  if(! utils::contained(id, _sendData)) {
    PtrCouplingData ptrCplData (new CouplingData(& (data->values()), initialize));
    DataMap::value_type pair = std::make_pair (id, ptrCplData);
    _sendData.insert(pair);
  }
  else {
    preciceError("addDataToSend()", "Data \"" << data->getName()
        << "\" of mesh \"" << data->mesh()->getName() << "\" cannot be "
        << "added twice for sending!");
  }
}

void BaseCouplingScheme:: addDataToReceive
(
   mesh::PtrData data,
   bool          initialize)
{
   int id = data->getID();
   if(! utils::contained(id, _receiveData)) {
      PtrCouplingData ptrCplData (new CouplingData(& (data->values()), initialize));
      DataMap::value_type pair = std::make_pair (id, ptrCplData);
      _receiveData.insert(pair);
   }
   else {
      preciceError("addDataToReceive()", "Data \"" << data->getName()
                     << "\" of mesh \"" << data->mesh()->getName() << "\" cannot be "
                     << "added twice for receiving!");
   }
}

void BaseCouplingScheme:: setCheckointTimestepInterval
(
  int timestepInterval)
{
  _checkpointTimestepInterval = timestepInterval;
}

void BaseCouplingScheme:: sendState
(
  com::PtrCommunication communication,
  int                   rankReceiver)
{
  preciceTrace("sendState()");
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
  foreach(const std::string& action, _actions){
    communication->send(action, rankReceiver);
  }
}

void BaseCouplingScheme:: receiveState
(
  com::PtrCommunication communication,
  int                   rankSender)
{
  preciceTrace("receiveState()");
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
  for (int i=0; i < actionsSize; i++){
    std::string action;
    communication->receive(action, rankSender);
    _actions.insert(action);
  }
}

std::vector<int> BaseCouplingScheme:: sendData
(
  com::PtrCommunication communication)
{
  preciceTrace("sendData()");
  assertion(communication.get() != NULL);
  assertion(communication->isConnected());

  std::vector<int> sentDataIDs;
  foreach (DataMap::value_type& pair, _sendData){
    int size = pair.second->values->size ();
    if (size > 0){
      communication->send(tarch::la::raw(*pair.second->values), size, 0);
    }
    sentDataIDs.push_back(pair.first);
  }
  preciceDebug("Number of sent data sets = " << sentDataIDs.size());
  return sentDataIDs;
}

std::vector<int> BaseCouplingScheme:: receiveData
(
  com::PtrCommunication communication)
{
  preciceTrace("receiveData()");
  assertion(communication.get() != NULL);
  assertion(communication->isConnected());

  std::vector<int> receivedDataIDs;
  foreach(DataMap::value_type & pair, _receiveData){
    int size = pair.second->values->size ();
    if (size > 0){
      communication->receive(tarch::la::raw(*pair.second->values), size, 0);
    }
    receivedDataIDs.push_back(pair.first);
  }
  preciceDebug("Number of received data sets = " << receivedDataIDs.size());
  return receivedDataIDs;
}

CouplingData* BaseCouplingScheme:: getSendData
(
  int dataID)
{
  preciceTrace1("getSendData()", dataID);
  DataMap::iterator iter = _sendData.find(dataID);
  if(iter != _sendData.end()){
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
  if(iter != _receiveData.end()){
    return  &(*(iter->second));
  }
  return NULL;
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

double BaseCouplingScheme:: getComputedTimestepPart() const
{
  return _computedTimestepPart;
}

void BaseCouplingScheme:: setComputedTimestepPart
(
  double computedTimestepPart)
{
  _computedTimestepPart = computedTimestepPart;
}

double BaseCouplingScheme:: getMaxTime() const
{
  return _maxTime;
}

int BaseCouplingScheme:: getMaxTimesteps() const
{
  return _maxTimesteps;
}

bool BaseCouplingScheme:: isInitialized() const
{
  return _isInitialized;
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
  double eps = std::pow(10.0, -1 * _validDigits);
  bool valid = tarch::la::greaterEquals(getThisTimestepRemainder(), 0.0, eps);
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
  double eps = std::pow(10.0, -1 * _validDigits);
  double remainder = getThisTimestepRemainder() - lastSolverTimestepLength;
  return not tarch::la::greater(remainder, 0.0, eps);
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

//double BaseCouplingScheme:: getTimestepLength() const
//{
//  if (tarch::la::equals(_timestepLength, UNDEFINED_TIMESTEP_LENGTH)){
//    return _computedTimestepPart;
//  }
//  return _timestepLength;
//}

double BaseCouplingScheme:: getThisTimestepRemainder() const
{
  preciceTrace("getTimestepRemainder()");
  double remainder = 0.0;
  if (not tarch::la::equals(_timestepLength, UNDEFINED_TIMESTEP_LENGTH)){
    remainder = _timestepLength - _computedTimestepPart;
  }
  preciceDebug("return " << remainder);
  return remainder;
//  preciceTrace1("getTimestepRemainder()", computedTimestepLength);
//  if(tarch::la::equals(_timestepLength, UNDEFINED_TIMESTEP_LENGTH)){
//    preciceDebug("Returning 0.0, since no timestep length is defined");
//    return 0.0;
//  }
//  double computedTimestepPart =
//      _time + computedTimestepLength - (_timesteps * _timestepLength);
//  double remainder = _timestepLength - computedTimestepPart;
//  double eps = std::pow(10.0, -1 * _validDigits);
//  preciceCheck(tarch::la::greaterEquals(remainder, 0.0, eps),
//                 "getTimestepRemainder()",
//                 "Computed timestep length (" << computedTimestepLength <<
//                 ") is not allowed to be larger than the prescribed one ("
//                 << remainder << ")! With " << _validDigits << " valid digits.");
//  return remainder;
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
  double eps = std::pow(10.0, -1 * _validDigits);
  using namespace tarch::la;
  bool timeLeft = greater(_maxTime, _time, eps) || equals(_maxTime, UNDEFINED_TIME);
  bool timestepsLeft = (_maxTimesteps > _timesteps)
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
   foreach(const std::string & actionName, _actions) {
      os << actionName << " | ";
   }
   return os.str ();
}

void BaseCouplingScheme:: checkCompletenessRequiredActions ()
{
  preciceTrace("checkCompletenessRequiredActions()");
  if(not _actions.empty()){
    std::ostringstream stream;
    foreach(const std::string & action, _actions){
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

}} // namespace precice, cplscheme
