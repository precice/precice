// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "CompositionalCouplingScheme.hpp"
#include "Constants.hpp"
#include "utils/Globals.hpp"
#include <limits>

namespace precice {
namespace cplscheme {

tarch::logging::Log CompositionalCouplingScheme::
   _log("precice::cplscheme::CompositionalCouplingScheme");

CompositionalCouplingScheme:: CompositionalCouplingScheme()
:
  _couplingSchemes(),
  _activeSchemesBegin(_couplingSchemes.begin()),
  _activeSchemesEnd(_couplingSchemes.end())
{}

void CompositionalCouplingScheme:: addCouplingScheme
(
  PtrCouplingScheme couplingScheme )
{
  preciceTrace("addCouplingScheme()");
  _couplingSchemes.push_back(couplingScheme);
}

void CompositionalCouplingScheme:: initialize
(
  double startTime,
  int    startTimestep )
{
  preciceTrace2("initialize()", startTime, startTimestep);
  foreach (PtrCouplingScheme couplingScheme, _couplingSchemes){
    couplingScheme->initialize(startTime, startTimestep);
  }
  determineActiveCouplingSchemes();
}

bool CompositionalCouplingScheme:: isInitialized() const
{
  preciceTrace("isInitialized()");
  bool isInitialized = true;
  foreach (PtrCouplingScheme couplingScheme, _couplingSchemes){
    isInitialized &= couplingScheme->isInitialized();
  }
  preciceDebug("return " << isInitialized);
  return isInitialized;
}

void CompositionalCouplingScheme:: initializeData()
{
  preciceTrace("initializeData()");
  foreach (PtrCouplingScheme couplingScheme, _couplingSchemes){
    couplingScheme->initializeData();
  }
}

void CompositionalCouplingScheme:: addComputedTime
(
  double timeToAdd )
{
  preciceTrace1("addComputedTime()", timeToAdd);
  for (SchemesIt it = _activeSchemesBegin; it != _activeSchemesEnd; it++){
    (*it)->addComputedTime(timeToAdd);
  }
}

void CompositionalCouplingScheme:: advance()
{
  preciceTrace("advance()");
  do {
    for (SchemesIt it = _activeSchemesBegin; it != _activeSchemesEnd; it++){
      (*it)->advance();
    }
  } while(determineActiveCouplingSchemes());
}

void CompositionalCouplingScheme:: finalize()
{
  preciceTrace("finalize()");
  foreach (PtrCouplingScheme couplingScheme, _couplingSchemes){
    couplingScheme->finalize();
  }
}

//bool CompositionalCouplingScheme:: isDataUsed
//(
//  int dataID )
//{
//  preciceTrace1("isDataUsed()", dataID);
//  bool isUsed = false;
//  foreach (PtrCouplingScheme couplingScheme, _couplingSchemes){
//    isUsed |= couplingScheme->isDataUsed(dataID);
//  }
//  return isUsed;
//}
//
//bool CompositionalCouplingScheme:: isDataUsed()
//{
//  preciceTrace("isDataUsed()");
//  bool isUsed = false;
//  foreach (PtrCouplingScheme couplingScheme, _couplingSchemes){
//    isUsed |= couplingScheme->isDataUsed();
//  }
//  return isUsed;
//}

std::vector<std::string> CompositionalCouplingScheme:: getCouplingPartners() const
{
  preciceTrace("getCouplingPartners()");
  std::vector<std::string> partners;
  std::vector<std::string> subpartners;
  foreach (PtrCouplingScheme couplingScheme, _couplingSchemes){
    subpartners = couplingScheme->getCouplingPartners();
    partners.insert(partners.end(), subpartners.begin(), subpartners.end());
  }
  return partners;
}

bool CompositionalCouplingScheme:: willDataBeExchanged
(
  double lastSolverTimestepLength ) const
{
  preciceTrace1("willDataBeExchanged()", lastSolverTimestepLength);
  bool willBeExchanged = false;
  for (SchemesIt it = _activeSchemesBegin; it != _activeSchemesEnd; it++){
    willBeExchanged |= (*it)->willDataBeExchanged(lastSolverTimestepLength);
  }
  preciceDebug("return " << willBeExchanged);
  return willBeExchanged;
}

bool CompositionalCouplingScheme:: hasDataBeenExchanged() const
{
  preciceTrace("hasDataBeenExchanged()");
  bool hasBeenExchanged = false;
  // Question: Does it suffice to only check the active ones?
  for (SchemesIt it = _activeSchemesBegin; it != _activeSchemesEnd; it++){
    hasBeenExchanged |= (*it)->hasDataBeenExchanged();
  }
  preciceDebug("return " << hasBeenExchanged);
  return hasBeenExchanged;
}

double CompositionalCouplingScheme:: getTime() const
{
  preciceTrace("getTime()");
  double time = std::numeric_limits<double>::max();
  foreach (PtrCouplingScheme couplingScheme, _couplingSchemes){
    if (couplingScheme->getTime() < time){
      time = couplingScheme->getTime();
    }
  }
  preciceDebug("return " << time);
  return time;
}

int CompositionalCouplingScheme:: getTimesteps() const
{
  preciceTrace("getTimesteps()");
  int timesteps = std::numeric_limits<int>::max();
  foreach (PtrCouplingScheme couplingScheme, _couplingSchemes){
    if (couplingScheme->getTimesteps() < timesteps){
      timesteps = couplingScheme->getTimesteps();
    }
  }
  preciceDebug("return " << timesteps);
  return timesteps;
}

double CompositionalCouplingScheme:: getMaxTime() const
{
  preciceTrace("getMaxTime()");
  double maxTime = 0.0;
  foreach (PtrCouplingScheme couplingScheme, _couplingSchemes){
    if (couplingScheme->getMaxTime() > maxTime){
      maxTime = couplingScheme->getMaxTime();
    }
  }
  preciceDebug("return " << maxTime);
  return maxTime;
}

int CompositionalCouplingScheme:: getMaxTimesteps() const
{
  preciceTrace("getMaxTimesteps()");
  int maxTimesteps = 0;
  foreach (PtrCouplingScheme couplingScheme, _couplingSchemes){
    if (couplingScheme->getMaxTimesteps() > maxTimesteps){
      maxTimesteps = couplingScheme->getMaxTimesteps();
    }
  }
  preciceDebug("return " << maxTimesteps);
  return maxTimesteps;
}

//int CompositionalCouplingScheme:: getSubIteration() const
//{
//  preciceTrace("getSubIteration()");
//  int subiteration = false;
//  foreach (PtrCouplingScheme couplingScheme, _couplingSchemes){
//    isSubiterating |= couplingScheme->getSubIteration();
//  }
//  return isSubiterating;
//}

bool CompositionalCouplingScheme:: hasTimestepLength() const
{
  preciceTrace("hasTimestepLength()");
  bool hasIt = false;
  foreach (PtrCouplingScheme couplingScheme, _couplingSchemes){
    hasIt |= couplingScheme->hasTimestepLength();
  }
  preciceDebug("return " << hasIt);
  return hasIt;
}

double CompositionalCouplingScheme:: getTimestepLength() const
{
  preciceTrace("getTimestepLength()");
  double timestepLength = std::numeric_limits<double>::max();
  foreach (PtrCouplingScheme couplingScheme, _couplingSchemes){
    if (couplingScheme->getTimestepLength() < timestepLength){
      timestepLength = couplingScheme->getTimestepLength();
    }
  }
  preciceDebug("return " << timestepLength);
  return timestepLength;
}

double CompositionalCouplingScheme:: getThisTimestepRemainder() const
{
  preciceTrace("getThisTimestepRemainder()");
  double maxRemainder = 0.0;
  foreach (PtrCouplingScheme couplingScheme, _couplingSchemes){
    if (couplingScheme->getThisTimestepRemainder() > maxRemainder){
      maxRemainder = couplingScheme->getThisTimestepRemainder();
    }
  }
  preciceDebug("return " << maxRemainder);
  return maxRemainder;
}

double CompositionalCouplingScheme:: getComputedTimestepPart() const
{
  preciceTrace("getComputedTimestepPart()");
  double timestepPart = std::numeric_limits<double>::max();
  foreach (PtrCouplingScheme couplingScheme, _couplingSchemes){
    if (couplingScheme->getComputedTimestepPart() < timestepPart){
      timestepPart = couplingScheme->getComputedTimestepPart();
    }
  }
  preciceDebug("return " << timestepPart);
  return timestepPart;
}

double CompositionalCouplingScheme:: getNextTimestepMaxLength() const
{
  preciceTrace("getNextTimestepMaxLength()");
  double maxLength = std::numeric_limits<double>::max();
  foreach (PtrCouplingScheme couplingScheme, _couplingSchemes){
    if (couplingScheme->getNextTimestepMaxLength() < maxLength){
      maxLength = couplingScheme->getNextTimestepMaxLength();
    }
  }
  preciceDebug("return " << maxLength);
  return maxLength;
}

bool CompositionalCouplingScheme:: isCouplingOngoing() const
{
  preciceTrace("isCouplingOngoing()");
  bool isOngoing = false;
  foreach (PtrCouplingScheme couplingScheme, _couplingSchemes){
    isOngoing |= couplingScheme->isCouplingOngoing();
  }
  preciceDebug("return " << isOngoing);
  return isOngoing;
}

bool CompositionalCouplingScheme:: isCouplingTimestepComplete() const
{
  preciceTrace("isCouplingTimestepComplete()");
  bool isComplete = true;
  foreach (PtrCouplingScheme couplingScheme, _couplingSchemes){
    isComplete &= couplingScheme->isCouplingTimestepComplete();
  }
  preciceDebug("return " << isComplete);
  return isComplete;
}

bool CompositionalCouplingScheme:: isActionRequired
(
  const std::string& actionName) const
{
  preciceTrace1("isActionRequired()", actionName);
  bool isRequired = false;
  foreach (PtrCouplingScheme couplingScheme, _couplingSchemes){
    isRequired |= couplingScheme->isActionRequired(actionName);
  }
  preciceDebug("return " << isRequired);
  return isRequired;
}

void CompositionalCouplingScheme:: performedAction
(
  const std::string& actionName)
{
  preciceTrace1("performedAction()", actionName);
  foreach (PtrCouplingScheme couplingScheme, _couplingSchemes){
    couplingScheme->performedAction(actionName);
  }
}

int CompositionalCouplingScheme:: getCheckpointTimestepInterval() const
{
  preciceTrace("getCheckpointTimestepInterval()");
  int interval = std::numeric_limits<int>::max();
  foreach (PtrCouplingScheme couplingScheme, _couplingSchemes){
    if (couplingScheme->getCheckpointTimestepInterval() < interval){
      interval = couplingScheme->getCheckpointTimestepInterval();
    }
  }
  preciceDebug("return " << interval);
  return interval;
}

void CompositionalCouplingScheme:: requireAction
(
  const std::string& actionName )
{
  preciceTrace1("requireAction()", actionName);
  foreach (PtrCouplingScheme couplingScheme, _couplingSchemes){
    couplingScheme->requireAction(actionName);
  }
}

std::string CompositionalCouplingScheme:: printCouplingState() const
{
  std::string state;
  std::vector<std::string> partners;
  foreach (PtrCouplingScheme couplingScheme, _couplingSchemes){
    if (not state.empty()){
      state += "\n";
    }
    partners = couplingScheme->getCouplingPartners();
    assertion1(partners.size() == 1, partners.size());
    state += "Coupling to ";
    state += partners[0];
    state += ":\n";
    state += couplingScheme->printCouplingState();
  }
  return state;
}

void CompositionalCouplingScheme:: exportState
(
  const std::string& filenamePrefix ) const
{
  preciceTrace("exportState()");
  int enumerator = 0;
  foreach (PtrCouplingScheme couplingScheme, _couplingSchemes){
    std::ostringstream stream;
    stream << filenamePrefix << "_" << enumerator;
    couplingScheme->exportState(stream.str());
    enumerator++;
  }
}

void CompositionalCouplingScheme:: importState
(
  const std::string& filenamePrefix )
{
  preciceTrace("importState()");
  int enumerator = 0;
  foreach (PtrCouplingScheme couplingScheme, _couplingSchemes){
    std::ostringstream stream;
    stream << filenamePrefix << "_" << enumerator;
    couplingScheme->importState(stream.str());
    enumerator++;
  }
}

void CompositionalCouplingScheme:: sendState
(
  com::PtrCommunication communication,
  int                   rankReceiver )
{
  preciceTrace("sendState()");
  foreach (PtrCouplingScheme couplingScheme, _couplingSchemes){
    couplingScheme->sendState(communication, rankReceiver);
  }
}

void CompositionalCouplingScheme:: receiveState
(
  com::PtrCommunication communication,
  int                   rankSender )
{
  preciceTrace("receiveState()");
  foreach (PtrCouplingScheme couplingScheme, _couplingSchemes){
    couplingScheme->receiveState(communication, rankSender);
  }
}

bool CompositionalCouplingScheme:: determineActiveCouplingSchemes()
{
  bool newActiveSchemes = false;
  std::string writeCheckpoint = constants::actionWriteIterationCheckpoint();
  std::string readCheckpoint = constants::actionReadIterationCheckpoint();
  if (_activeSchemesBegin == _activeSchemesEnd){
    // First call after initialization of all coupling schemes. All coupling
    // schemes are set active up to (but not including) the first explicit
    // scheme after an implicit scheme.
    _activeSchemesBegin = _couplingSchemes.begin();
    _activeSchemesEnd = _couplingSchemes.begin();
    advanceActiveCouplingSchemes();
    newActiveSchemes = true;
  }
  else {
    // Redetermine active schemes. First, all preceding explicit schemes are
    // removed. Then, all remaining implicit schemes are checked for the
    // convergence of iterations (this is given when an iteration checkpoint
    // should be created). If all are converged, a next set of active schemes
    // is determined.

    // Remove preceding explicit schemes
    while (_activeSchemesBegin != _activeSchemesEnd){
      bool explicitScheme = true;
      explicitScheme &= not (*_activeSchemesBegin)->isActionRequired(writeCheckpoint);
      explicitScheme &= not (*_activeSchemesBegin)->isActionRequired(readCheckpoint);
      if (explicitScheme) _activeSchemesBegin++;
    }

    // Check implicit schemes for convergence and remove if converged
    bool converged = true;
    for (SchemesIt it=_activeSchemesBegin; it != _activeSchemesEnd; it++){
      if ((*it)->isActionRequired(readCheckpoint)){
        converged = false;
        break;
      }
    }
    if (converged) _activeSchemesBegin = _activeSchemesEnd;

    // Determine next set of active schemes if current is empty
    if (_activeSchemesBegin == _activeSchemesEnd){
      if (_activeSchemesBegin == _couplingSchemes.end()){
        // All coupling schemes are through
        _activeSchemesBegin = _couplingSchemes.begin();
        _activeSchemesEnd = _couplingSchemes.begin();
        advanceActiveCouplingSchemes();
        // newActiveSchemes stays false, since the current it/dt is complete
      }
      else {
        advanceActiveCouplingSchemes();
        newActiveSchemes = true;
      }
    }
  }
  return newActiveSchemes;
}

void CompositionalCouplingScheme:: advanceActiveCouplingSchemes()
{
  std::string writeCheckpoint = constants::actionWriteIterationCheckpoint();
  bool iterating = false;
  while (_activeSchemesEnd != _couplingSchemes.end()){
    if ((*_activeSchemesEnd)->isActionRequired(writeCheckpoint)){
      iterating = true;
    }
    if (iterating && (not (*_activeSchemesEnd)->isActionRequired(writeCheckpoint))){
      break;
    }
    _activeSchemesEnd++;
  }
  assertion(_activeSchemesBegin != _activeSchemesEnd);
}

}} // namespace precice, cplscheme
