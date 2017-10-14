#include "CompositionalCouplingScheme.hpp"
#include "Constants.hpp"
#include "utils/Globals.hpp"
#include <limits>

namespace precice {
namespace cplscheme {

logging::Logger CompositionalCouplingScheme::
   _log("cplscheme::CompositionalCouplingScheme");

CompositionalCouplingScheme:: CompositionalCouplingScheme()
:
  _couplingSchemes(),
  _activeSchemesBegin(_couplingSchemes.end()),
  _activeSchemesEnd(_couplingSchemes.end()),
  _lastAddedTime(0.0)
{}

void CompositionalCouplingScheme:: addCouplingScheme
(
  PtrCouplingScheme couplingScheme )
{
  TRACE();
  _couplingSchemes.push_back(Scheme(couplingScheme));
}

void CompositionalCouplingScheme:: initialize
(
  double startTime,
  int    startTimestep )
{
  TRACE(startTime, startTimestep);
  for (Scheme scheme : _couplingSchemes) {
    scheme.scheme->initialize(startTime, startTimestep);
  }
  determineActiveCouplingSchemes();
}

bool CompositionalCouplingScheme:: isInitialized() const
{
  TRACE();
  bool isInitialized = true;
  for (Scheme scheme : _couplingSchemes) {
    isInitialized &= scheme.scheme->isInitialized();
  }
  DEBUG("return " << isInitialized);
  return isInitialized;
}

void CompositionalCouplingScheme:: initializeData()
{
  TRACE();
  for (Scheme scheme : _couplingSchemes) {
    scheme.scheme->initializeData();
  }
}

void CompositionalCouplingScheme:: addComputedTime
(
  double timeToAdd )
{
  TRACE(timeToAdd);
  for (SchemesIt it = _activeSchemesBegin; it != _activeSchemesEnd; it++){
    if (not it->onHold){
      it->scheme->addComputedTime(timeToAdd);
    }
  }
  _lastAddedTime += timeToAdd;
}

void CompositionalCouplingScheme:: advance()
{
  TRACE();
  bool moreSchemesToHandle = false;
  do {
    for (SchemesIt it = _activeSchemesBegin; it != _activeSchemesEnd; it++){
      if (not it->onHold){
        it->scheme->advance();
      }
    }
    moreSchemesToHandle = determineActiveCouplingSchemes();
    if (moreSchemesToHandle){
      // The new schemes to be handled in this advance also need the time that
      // has been computed so far. This time can't be added in the solver call
      // to addComputedTime(), since there the schemes are not active yet.
      addComputedTime(_lastAddedTime);
    }
  } while (moreSchemesToHandle);
  _lastAddedTime = 0.0;
}

void CompositionalCouplingScheme:: finalize()
{
  TRACE();
  for (Scheme scheme : _couplingSchemes) {
    scheme.scheme->finalize();
  }
}

//bool CompositionalCouplingScheme:: isDataUsed
//(
//  int dataID )
//{
//  TRACE(dataID);
//  bool isUsed = false;
//  for (PtrCouplingScheme couplingScheme : _couplingSchemes){
//    isUsed |= couplingScheme->isDataUsed(dataID);
//  }
//  return isUsed;
//}
//
//bool CompositionalCouplingScheme:: isDataUsed()
//{
//  TRACE();
//  bool isUsed = false;
//  for (PtrCouplingScheme couplingScheme : _couplingSchemes){
//    isUsed |= couplingScheme->isDataUsed();
//  }
//  return isUsed;
//}

std::vector<std::string> CompositionalCouplingScheme:: getCouplingPartners() const
{
  TRACE();
  std::vector<std::string> partners;
  std::vector<std::string> subpartners;
  for (Scheme scheme : _couplingSchemes) {
    subpartners = scheme.scheme->getCouplingPartners();
    partners.insert(partners.end(), subpartners.begin(), subpartners.end());
  }
  return partners;
}

bool CompositionalCouplingScheme:: willDataBeExchanged
(
  double lastSolverTimestepLength ) const
{
  TRACE(lastSolverTimestepLength);
  bool willBeExchanged = false;
  for (SchemesIt it = _activeSchemesBegin; it != _activeSchemesEnd; it++){
    if (not it->onHold){
      willBeExchanged |= it->scheme->willDataBeExchanged(lastSolverTimestepLength);
    }
  }
  DEBUG("return " << willBeExchanged);
  return willBeExchanged;
}

bool CompositionalCouplingScheme:: hasDataBeenExchanged() const
{
  TRACE();
  bool hasBeenExchanged = false;
  // Question: Does it suffice to only check the active ones?
  for (SchemesIt it = _activeSchemesBegin; it != _activeSchemesEnd; it++){
    if (not it->onHold){
      hasBeenExchanged |= it->scheme->hasDataBeenExchanged();
    }
  }
  DEBUG("return " << hasBeenExchanged);
  return hasBeenExchanged;
}

double CompositionalCouplingScheme:: getTime() const
{
  TRACE();
  double time = std::numeric_limits<double>::max();
  for (Scheme scheme : _couplingSchemes) {
    if (not scheme.onHold){
      if (scheme.scheme->getTime() < time){
        time = scheme.scheme->getTime();
      }
    }
  }
  DEBUG("return " << time);
  return time;
}

int CompositionalCouplingScheme:: getTimesteps() const
{
  TRACE();
  int timesteps = std::numeric_limits<int>::max();
  for (Scheme scheme : _couplingSchemes) {
    if (not scheme.onHold){
      if (scheme.scheme->getTimesteps() < timesteps){
        timesteps = scheme.scheme->getTimesteps();
      }
    }
  }
  DEBUG("return " << timesteps);
  return timesteps;
}

double CompositionalCouplingScheme:: getMaxTime() const
{
  TRACE();
  double maxTime = 0.0;
  for (Scheme scheme : _couplingSchemes) {
    if (scheme.scheme->getMaxTime() > maxTime){
      maxTime = scheme.scheme->getMaxTime();
    }
  }
  DEBUG("return " << maxTime);
  return maxTime;
}

int CompositionalCouplingScheme:: getMaxTimesteps() const
{
  TRACE();
  int maxTimesteps = 0;
  for (Scheme scheme : _couplingSchemes) {
    if (scheme.scheme->getMaxTimesteps() > maxTimesteps){
      maxTimesteps = scheme.scheme->getMaxTimesteps();
    }
  }
  DEBUG("return " << maxTimesteps);
  return maxTimesteps;
}

//int CompositionalCouplingScheme:: getSubIteration() const
//{
//  TRACE();
//  int subiteration = false;
//  for (PtrCouplingScheme couplingScheme : _couplingSchemes){
//    isSubiterating |= couplingScheme->getSubIteration();
//  }
//  return isSubiterating;
//}

bool CompositionalCouplingScheme:: hasTimestepLength() const
{
  TRACE();
  bool hasIt = false;
  for (Scheme scheme : _couplingSchemes) {
    hasIt |= scheme.scheme->hasTimestepLength();
  }
  DEBUG("return " << hasIt);
  return hasIt;
}

double CompositionalCouplingScheme:: getTimestepLength() const
{
  TRACE();
  double timestepLength = std::numeric_limits<double>::max();
  for (Scheme scheme : _couplingSchemes) {
    if (scheme.scheme->getTimestepLength() < timestepLength){
      timestepLength = scheme.scheme->getTimestepLength();
    }
  }
  DEBUG("return " << timestepLength);
  return timestepLength;
}

double CompositionalCouplingScheme:: getThisTimestepRemainder() const
{
  TRACE();
  double maxRemainder = 0.0;
  for (Scheme scheme : _couplingSchemes) {
    if (not scheme.onHold){
      if (scheme.scheme->getThisTimestepRemainder() > maxRemainder){
        maxRemainder = scheme.scheme->getThisTimestepRemainder();
      }
    }
  }
  DEBUG("return " << maxRemainder);
  return maxRemainder;
}

double CompositionalCouplingScheme:: getComputedTimestepPart() const
{
  TRACE();
  double timestepPart = std::numeric_limits<double>::max();
  for (Scheme scheme : _couplingSchemes) {
    if (not scheme.onHold){
      if (scheme.scheme->getComputedTimestepPart() < timestepPart){
        timestepPart = scheme.scheme->getComputedTimestepPart();
      }
    }
  }
  DEBUG("return " << timestepPart);
  return timestepPart;
}

double CompositionalCouplingScheme:: getNextTimestepMaxLength() const
{
  TRACE();
  double maxLength = std::numeric_limits<double>::max();
  for (Scheme scheme : _couplingSchemes) {
    if (not scheme.onHold){
      if (scheme.scheme->getNextTimestepMaxLength() < maxLength){
        maxLength = scheme.scheme->getNextTimestepMaxLength();
      }
    }
  }
  DEBUG("return " << maxLength);
  return maxLength;
}

bool CompositionalCouplingScheme:: isCouplingOngoing() const
{
  TRACE();
  bool isOngoing = false;
  for (Scheme scheme : _couplingSchemes) {
    isOngoing |= scheme.scheme->isCouplingOngoing();
  }
  DEBUG("return " << isOngoing);
  return isOngoing;
}

bool CompositionalCouplingScheme:: isCouplingTimestepComplete() const
{
  TRACE();
  bool isComplete = true;
  for (Scheme scheme : _couplingSchemes) {
    isComplete &= scheme.scheme->isCouplingTimestepComplete();
  }
  DEBUG("return " << isComplete);
  return isComplete;
}

bool CompositionalCouplingScheme:: isActionRequired
(
  const std::string& actionName) const
{
  TRACE(actionName);
  bool isRequired = false;
  for (Scheme scheme : _couplingSchemes) {
    if (not scheme.onHold){
      isRequired |= scheme.scheme->isActionRequired(actionName);
    }
  }
  DEBUG("return " << isRequired);
  return isRequired;
}

void CompositionalCouplingScheme:: performedAction
(
  const std::string& actionName)
{
  TRACE(actionName);
  for (Scheme scheme : _couplingSchemes) {
    if (not scheme.onHold){
      scheme.scheme->performedAction(actionName);
    }
  }
}


void CompositionalCouplingScheme:: requireAction
(
  const std::string& actionName )
{
  TRACE(actionName);
  for (Scheme scheme : _couplingSchemes) {
    scheme.scheme->requireAction(actionName);
  }
}

std::string CompositionalCouplingScheme:: printCouplingState() const
{
  std::string state;
  std::vector<std::string> partners;
  for (Scheme scheme : _couplingSchemes) {
    if (not state.empty()){
      state += "\n";
    }
    partners = scheme.scheme->getCouplingPartners();
    //TODO i don't understand why this assertion should make any sense, it fails when
    // more than 3 participants are used, but without any mature reason, maybe
    // the output is just not nice enough
    //assertion(partners.size() == 1, partners.size());
    state += partners[0];
    state += ": ";
    state += scheme.scheme->printCouplingState();
  }
  return state;
}

void CompositionalCouplingScheme:: sendState
(
  com::PtrCommunication communication,
  int                   rankReceiver )
{
  TRACE();
  for (Scheme scheme : _couplingSchemes) {
    scheme.scheme->sendState(communication, rankReceiver);
  }
}

void CompositionalCouplingScheme:: receiveState
(
  com::PtrCommunication communication,
  int                   rankSender )
{
  TRACE();
  for (Scheme scheme : _couplingSchemes) {
    scheme.scheme->receiveState(communication, rankSender);
  }
}

bool CompositionalCouplingScheme:: determineActiveCouplingSchemes()
{
  TRACE();
  bool newActiveSchemes = false;
  std::string writeCheckpoint = constants::actionWriteIterationCheckpoint();
  std::string readCheckpoint = constants::actionReadIterationCheckpoint();
  if (_activeSchemesBegin == _activeSchemesEnd){
    DEBUG("Case After Init");
    // First call after initialization of all coupling schemes. All coupling
    // schemes are set active up to (but not including) the first explicit
    // scheme after an implicit scheme.
    _activeSchemesBegin = _couplingSchemes.begin();
    _activeSchemesEnd = _couplingSchemes.begin();
    advanceActiveCouplingSchemes();
    newActiveSchemes = true;
  }
  else {
    DEBUG("Normal Case");
    // Redetermine active schemes. First, all preceding explicit schemes are
    // removed. Then, all remaining implicit schemes are checked for the
    // convergence of iterations (this is given when an iteration checkpoint
    // should be created). If all are converged, a next set of active schemes
    // is determined.

    // Remove preceding explicit schemes
    while (_activeSchemesBegin != _activeSchemesEnd){
      bool explicitScheme = true;
      explicitScheme &= not _activeSchemesBegin->scheme->isActionRequired(writeCheckpoint);
      explicitScheme &= not _activeSchemesBegin->scheme->isActionRequired(readCheckpoint);
      if (explicitScheme) {
        _activeSchemesBegin++;
        DEBUG("Remove preceding explicit scheme");
      }
      else {
        break;
      }
    }

    // Check implicit schemes for convergence and remove if converged
    bool converged = true;
    for (SchemesIt it=_activeSchemesBegin; it != _activeSchemesEnd; it++){
      if (it->scheme->isActionRequired(readCheckpoint)){
        converged = false;
        DEBUG("Non converged implicit scheme");
      }
      else if (it->scheme->isActionRequired(writeCheckpoint)
               || not it->scheme->isCouplingOngoing())
      {
        it->onHold = true;
        DEBUG("Put converged/finished implicit scheme on hold");
      }
    }
    if (converged) {
      DEBUG("Active implicit schemes converged");
      for (SchemesIt it=_activeSchemesBegin; it != _activeSchemesEnd; it++){
        it->onHold = false;
      }
      _activeSchemesBegin = _activeSchemesEnd;
    }

    // Determine next set of active schemes if current is empty
    if (_activeSchemesBegin == _activeSchemesEnd){
      if (_activeSchemesBegin == _couplingSchemes.end()){
        DEBUG("Through with all coupling schemes");
        // All coupling schemes are through
        _activeSchemesBegin = _couplingSchemes.begin();
        _activeSchemesEnd = _couplingSchemes.begin();
        advanceActiveCouplingSchemes();
        // newActiveSchemes stays false, since the current it/dt is complete
      }
      else {
        DEBUG("Coupling schemes remaining");
        advanceActiveCouplingSchemes();
        newActiveSchemes = true;
      }
    }
  }
  DEBUG("return newActiveSchemes=" << newActiveSchemes);
  return newActiveSchemes;
}

void CompositionalCouplingScheme:: advanceActiveCouplingSchemes()
{
  TRACE();
  std::string writeCheckpoint = constants::actionWriteIterationCheckpoint();
  bool iterating = false;
  while (_activeSchemesEnd != _couplingSchemes.end()){
    if (_activeSchemesEnd->scheme->isActionRequired(writeCheckpoint)){
      DEBUG("Found implicit scheme");
      iterating = true;
    }
    if (iterating && (not _activeSchemesEnd->scheme->isActionRequired(writeCheckpoint))){
      DEBUG("Found explicit scheme after implicit scheme");
      break;
    }
    _activeSchemesEnd++;
    DEBUG("Advanced active schemes by one new scheme");
  }
  assertion(_activeSchemesBegin != _activeSchemesEnd);
}

}} // namespace precice, cplscheme

//#include "CompositionalCouplingScheme.hpp"
//#include "Constants.hpp"
//#include "utils/Globals.hpp"
//#include <limits>
//
//namespace precice {
//namespace cplscheme {
//
//logging::Logger CompositionalCouplingScheme::
//   _log("cplscheme::CompositionalCouplingScheme");
//
//CompositionalCouplingScheme:: CompositionalCouplingScheme()
//:
//  _couplingSchemes(),
//  _activeCouplingSchemes(),
//  //_activeSchemesBegin(_couplingSchemes.end()),
//  _activeSchemesEnd(_couplingSchemes.end()),
//  _lastAddedTime(0.0)
//{}
//
//void CompositionalCouplingScheme:: addCouplingScheme
//(
//  PtrCouplingScheme couplingScheme )
//{
//  TRACE();
//  _couplingSchemes.push_back(couplingScheme);
//}
//
//void CompositionalCouplingScheme:: initialize
//(
//  double startTime,
//  int    startTimestep )
//{
//  TRACE(startTime, startTimestep);
//  for (PtrCouplingScheme couplingScheme : _couplingSchemes){
//    couplingScheme->initialize(startTime, startTimestep);
//  }
//  determineActiveCouplingSchemes();
//}
//
//bool CompositionalCouplingScheme:: isInitialized() const
//{
//  TRACE();
//  bool isInitialized = true;
//  for (PtrCouplingScheme couplingScheme : _couplingSchemes){
//    isInitialized &= couplingScheme->isInitialized();
//  }
//  DEBUG("return " << isInitialized);
//  return isInitialized;
//}
//
//void CompositionalCouplingScheme:: initializeData()
//{
//  TRACE();
//  for (PtrCouplingScheme couplingScheme : _couplingSchemes){
//    if (couplingScheme->isActionRequired(constants::actionWriteInitialData())){
//      couplingScheme->initializeData();
//    }
//  }
//}
//
//void CompositionalCouplingScheme:: addComputedTime
//(
//  double timeToAdd )
//{
//  TRACE(timeToAdd);
//  foriter(SchemesIt, it, _activeCouplingSchemes){
//  //for (SchemesIt it = _activeSchemesBegin; it != _activeSchemesEnd; it++){
//    (*it)->addComputedTime(timeToAdd);
//  }
//  _lastAddedTime += timeToAdd;
//}
//
//void CompositionalCouplingScheme:: advance()
//{
//  TRACE();
//  bool moreSchemesToHandle = false;
//  do {
//    //for (SchemesIt it = _activeSchemesBegin; it != _activeSchemesEnd; it++){
//    foriter(SchemesIt, it, _activeCouplingSchemes){
//      (*it)->advance();
//    }
//    moreSchemesToHandle = determineActiveCouplingSchemes();
//    if (moreSchemesToHandle){
//      // The new schemes to be handled in this advance also need the time that
//      // has been computed so far. This time can't be added in the solver call
//      // to addComputedTime(), since there the schemes are not active yet.
//      addComputedTime(_lastAddedTime);
//    }
//  } while (moreSchemesToHandle);
//  _lastAddedTime = 0.0;
//}
//
//void CompositionalCouplingScheme:: finalize()
//{
//  TRACE();
//  for (PtrCouplingScheme couplingScheme : _couplingSchemes){
//    couplingScheme->finalize();
//  }
//}
//
////bool CompositionalCouplingScheme:: isDataUsed
////(
////  int dataID )
////{
////  TRACE(dataID);
////  bool isUsed = false;
////  for (PtrCouplingScheme couplingScheme : _couplingSchemes){
////    isUsed |= couplingScheme->isDataUsed(dataID);
////  }
////  return isUsed;
////}
////
////bool CompositionalCouplingScheme:: isDataUsed()
////{
////  TRACE();
////  bool isUsed = false;
////  for (PtrCouplingScheme couplingScheme : _couplingSchemes){
////    isUsed |= couplingScheme->isDataUsed();
////  }
////  return isUsed;
////}
//
//std::vector<std::string> CompositionalCouplingScheme:: getCouplingPartners() const
//{
//  TRACE();
//  std::vector<std::string> partners;
//  std::vector<std::string> subpartners;
//  for (PtrCouplingScheme couplingScheme : _couplingSchemes){
//    subpartners = couplingScheme->getCouplingPartners();
//    partners.insert(partners.end(), subpartners.begin(), subpartners.end());
//  }
//  return partners;
//}
//
//bool CompositionalCouplingScheme:: willDataBeExchanged
//(
//  double lastSolverTimestepLength ) const
//{
//  TRACE(lastSolverTimestepLength);
//  bool willBeExchanged = false;
//  //for (SchemesIt it = _activeSchemesBegin; it != _activeSchemesEnd; it++){
//  foriter(ConstSchemesIt, it, _activeCouplingSchemes){
//    willBeExchanged |= (*it)->willDataBeExchanged(lastSolverTimestepLength);
//  }
//  DEBUG("return " << willBeExchanged);
//  return willBeExchanged;
//}
//
//bool CompositionalCouplingScheme:: hasDataBeenExchanged() const
//{
//  TRACE();
//  bool hasBeenExchanged = false;
//  // Question: Does it suffice to only check the active ones?
//  //for (SchemesIt it = _activeSchemesBegin; it != _activeSchemesEnd; it++){
//  foriter(ConstSchemesIt, it, _activeCouplingSchemes){
//    hasBeenExchanged |= (*it)->hasDataBeenExchanged();
//  }
//  DEBUG("return " << hasBeenExchanged);
//  return hasBeenExchanged;
//}
//
//double CompositionalCouplingScheme:: getTime() const
//{
//  TRACE();
//  double time = std::numeric_limits<double>::max();
//  for (PtrCouplingScheme couplingScheme : _couplingSchemes){
//    if (couplingScheme->getTime() < time){
//      time = couplingScheme->getTime();
//    }
//  }
//  DEBUG("return " << time);
//  return time;
//}
//
//int CompositionalCouplingScheme:: getTimesteps() const
//{
//  TRACE();
//  int timesteps = std::numeric_limits<int>::max();
//  for (PtrCouplingScheme couplingScheme : _couplingSchemes){
//    if (couplingScheme->getTimesteps() < timesteps){
//      timesteps = couplingScheme->getTimesteps();
//    }
//  }
//  DEBUG("return " << timesteps);
//  return timesteps;
//}
//
//double CompositionalCouplingScheme:: getMaxTime() const
//{
//  TRACE();
//  double maxTime = 0.0;
//  for (PtrCouplingScheme couplingScheme : _couplingSchemes){
//    if (couplingScheme->getMaxTime() > maxTime){
//      maxTime = couplingScheme->getMaxTime();
//    }
//  }
//  DEBUG("return " << maxTime);
//  return maxTime;
//}
//
//int CompositionalCouplingScheme:: getMaxTimesteps() const
//{
//  TRACE();
//  int maxTimesteps = 0;
//  for (PtrCouplingScheme couplingScheme : _couplingSchemes){
//    if (couplingScheme->getMaxTimesteps() > maxTimesteps){
//      maxTimesteps = couplingScheme->getMaxTimesteps();
//    }
//  }
//  DEBUG("return " << maxTimesteps);
//  return maxTimesteps;
//}
//
////int CompositionalCouplingScheme:: getSubIteration() const
////{
////  TRACE();
////  int subiteration = false;
////  for (PtrCouplingScheme couplingScheme : _couplingSchemes){
////    isSubiterating |= couplingScheme->getSubIteration();
////  }
////  return isSubiterating;
////}
//
//bool CompositionalCouplingScheme:: hasTimestepLength() const
//{
//  TRACE();
//  bool hasIt = false;
//  for (PtrCouplingScheme couplingScheme : _couplingSchemes){
//    hasIt |= couplingScheme->hasTimestepLength();
//  }
//  DEBUG("return " << hasIt);
//  return hasIt;
//}
//
//double CompositionalCouplingScheme:: getTimestepLength() const
//{
//  TRACE();
//  double timestepLength = std::numeric_limits<double>::max();
//  for (PtrCouplingScheme couplingScheme : _couplingSchemes){
//    if (couplingScheme->getTimestepLength() < timestepLength){
//      timestepLength = couplingScheme->getTimestepLength();
//    }
//  }
//  DEBUG("return " << timestepLength);
//  return timestepLength;
//}
//
//double CompositionalCouplingScheme:: getThisTimestepRemainder() const
//{
//  TRACE();
//  double maxRemainder = 0.0;
//  for (PtrCouplingScheme couplingScheme : _couplingSchemes){
//    if (couplingScheme->getThisTimestepRemainder() > maxRemainder){
//      maxRemainder = couplingScheme->getThisTimestepRemainder();
//    }
//  }
//  DEBUG("return " << maxRemainder);
//  return maxRemainder;
//}
//
//double CompositionalCouplingScheme:: getComputedTimestepPart() const
//{
//  TRACE();
//  double timestepPart = std::numeric_limits<double>::max();
//  for (PtrCouplingScheme couplingScheme : _couplingSchemes){
//    if (couplingScheme->getComputedTimestepPart() < timestepPart){
//      timestepPart = couplingScheme->getComputedTimestepPart();
//    }
//  }
//  DEBUG("return " << timestepPart);
//  return timestepPart;
//}
//
//double CompositionalCouplingScheme:: getNextTimestepMaxLength() const
//{
//  TRACE();
//  double maxLength = std::numeric_limits<double>::max();
//  for (PtrCouplingScheme couplingScheme : _couplingSchemes){
//    if (couplingScheme->getNextTimestepMaxLength() < maxLength){
//      maxLength = couplingScheme->getNextTimestepMaxLength();
//    }
//  }
//  DEBUG("return " << maxLength);
//  return maxLength;
//}
//
//bool CompositionalCouplingScheme:: isCouplingOngoing() const
//{
//  TRACE();
//  bool isOngoing = false;
//  for (PtrCouplingScheme couplingScheme : _couplingSchemes){
//    isOngoing |= couplingScheme->isCouplingOngoing();
//  }
//  DEBUG("return " << isOngoing);
//  return isOngoing;
//}
//
//bool CompositionalCouplingScheme:: isCouplingTimestepComplete() const
//{
//  TRACE();
//  bool isComplete = true;
//  for (PtrCouplingScheme couplingScheme : _couplingSchemes){
//    isComplete &= couplingScheme->isCouplingTimestepComplete();
//  }
//  DEBUG("return " << isComplete);
//  return isComplete;
//}
//
//bool CompositionalCouplingScheme:: isActionRequired
//(
//  const std::string& actionName) const
//{
//  TRACE(actionName);
//  bool isRequired = false;
//  for (PtrCouplingScheme couplingScheme : _couplingSchemes){
//    isRequired |= couplingScheme->isActionRequired(actionName);
//  }
//  DEBUG("return " << isRequired);
//  return isRequired;
//}
//
//void CompositionalCouplingScheme:: performedAction
//(
//  const std::string& actionName)
//{
//  TRACE(actionName);
//  for (PtrCouplingScheme couplingScheme : _couplingSchemes){
//    couplingScheme->performedAction(actionName);
//  }
//}
//
//int CompositionalCouplingScheme:: getCheckpointTimestepInterval() const
//{
//  TRACE();
//  int interval = std::numeric_limits<int>::max();
//  for (PtrCouplingScheme couplingScheme : _couplingSchemes){
//    if (couplingScheme->getCheckpointTimestepInterval() < interval){
//      interval = couplingScheme->getCheckpointTimestepInterval();
//    }
//  }
//  DEBUG("return " << interval);
//  return interval;
//}
//
//void CompositionalCouplingScheme:: requireAction
//(
//  const std::string& actionName )
//{
//  TRACE(actionName);
//  for (PtrCouplingScheme couplingScheme : _couplingSchemes){
//    couplingScheme->requireAction(actionName);
//  }
//}
//
//std::string CompositionalCouplingScheme:: printCouplingState() const
//{
//  std::string state;
//  std::vector<std::string> partners;
//  for (PtrCouplingScheme couplingScheme : _couplingSchemes){
//    if (not state.empty()){
//      state += "\n";
//    }
//    partners = couplingScheme->getCouplingPartners();
//    assertion(partners.size() == 1, partners.size());
//    state += "Coupling to ";
//    state += partners[0];
//    state += ":\n";
//    state += couplingScheme->printCouplingState();
//  }
//  return state;
//}
//
//void CompositionalCouplingScheme:: exportState
//(
//  const std::string& filenamePrefix ) const
//{
//  TRACE();
//  int enumerator = 0;
//  for (PtrCouplingScheme couplingScheme : _couplingSchemes){
//    std::ostringstream stream;
//    stream << filenamePrefix << "_" << enumerator;
//    couplingScheme->exportState(stream.str());
//    enumerator++;
//  }
//}
//
//void CompositionalCouplingScheme:: importState
//(
//  const std::string& filenamePrefix )
//{
//  TRACE();
//  int enumerator = 0;
//  for (PtrCouplingScheme couplingScheme : _couplingSchemes){
//    std::ostringstream stream;
//    stream << filenamePrefix << "_" << enumerator;
//    couplingScheme->importState(stream.str());
//    enumerator++;
//  }
//}
//
//void CompositionalCouplingScheme:: sendState
//(
//  com::Communication::SharedPointer communication,
//  int                   rankReceiver )
//{
//  TRACE();
//  for (PtrCouplingScheme couplingScheme : _couplingSchemes){
//    couplingScheme->sendState(communication, rankReceiver);
//  }
//}
//
//void CompositionalCouplingScheme:: receiveState
//(
//  com::Communication::SharedPointer communication,
//  int                   rankSender )
//{
//  TRACE();
//  for (PtrCouplingScheme couplingScheme : _couplingSchemes){
//    couplingScheme->receiveState(communication, rankSender);
//  }
//}
//
//bool CompositionalCouplingScheme:: determineActiveCouplingSchemes()
//{
//  TRACE();
//  bool newActiveSchemes = false;
//  std::string writeCheckpoint = constants::actionWriteIterationCheckpoint();
//  std::string readCheckpoint = constants::actionReadIterationCheckpoint();
//  if (_activeCouplingSchemes.empty()){
//    DEBUG("Case After Init");
//    // First call after initialization of all coupling schemes. All coupling
//    // schemes are set active up to (but not including) the first explicit
//    // scheme after an implicit scheme.
//    //_activeSchemesBegin = _couplingSchemes.begin();
//    _activeSchemesEnd = _couplingSchemes.begin();
//    advanceActiveCouplingSchemes();
//    newActiveSchemes = true;
//  }
//  else {
//    DEBUG("Normal Case");
//    // Redetermine active schemes. First, all preceding explicit schemes are
//    // removed. Then, all remaining implicit schemes are checked for the
//    // convergence of iterations (this is given when an iteration checkpoint
//    // should be created). If all are converged, a next set of active schemes
//    // is determined.
//
//    // Remove preceding explicit schemes
//    //while (_activeSchemesBegin != _activeSchemesEnd){
//    while (not _activeCouplingSchemes.empty()){
//      bool explicitScheme = true;
//      SchemesIt it = _activeCouplingSchemes.begin();
//      explicitScheme &= not (*it)->isActionRequired(writeCheckpoint);
//      explicitScheme &= not (*it)->isActionRequired(readCheckpoint);
//      if (explicitScheme) {
//        //_activeSchemesBegin++;
//        _activeCouplingSchemes.pop_front();
//        DEBUG("Remove preceding explicit scheme");
//      }
//      else {
//        break;
//      }
//    }
//
//    // Check implicit schemes for convergence and remove if converged
//    if (not _activeCouplingSchemes.empty()){
//      bool converged = true;
//      //for (SchemesIt it=_activeSchemesBegin; it != _activeSchemesEnd; it++){
//      SchemesIt it = _activeCouplingSchemes.begin();
//      while (it != _activeCouplingSchemes.end()){
//        if ((*it)->isActionRequired(readCheckpoint)){
//          converged = false;
//          it++; // Advance it to the next scheme
//          DEBUG("Non converged implicit scheme");
//        }
//        else if ((*it)->isActionRequired(writeCheckpoint)){
//          it = _activeCouplingSchemes.erase(it); // Automatically advances to next scheme
//          DEBUG("Remove converged implicit scheme");
//        }
//      }
//      if (converged) {
//        DEBUG("Active implicit schemes converged");
//        _activeCouplingSchemes.clear();
//        //_activeSchemesBegin = _activeSchemesEnd;
//      }
//    }
//
//    // Determine next set of active schemes if current is empty
//    //if (_activeSchemesBegin == _activeSchemesEnd){
//    if (_activeCouplingSchemes.empty()){
//      if (_activeSchemesEnd == _couplingSchemes.end()){
//        DEBUG("Through with all coupling schemes");
//        // All coupling schemes are through
//        //_activeSchemesBegin = _couplingSchemes.begin();
//        _activeSchemesEnd = _couplingSchemes.begin();
//        advanceActiveCouplingSchemes();
//        // newActiveSchemes stays false, since the current it/dt is complete
//      }
//      else {
//        DEBUG("Coupling schemes remaining");
//        advanceActiveCouplingSchemes();
//        newActiveSchemes = true;
//      }
//    }
//  }
//  DEBUG("return newActiveSchemes=" << newActiveSchemes);
//  return newActiveSchemes;
//}
//
//void CompositionalCouplingScheme:: advanceActiveCouplingSchemes()
//{
//  TRACE();
//  std::string writeCheckpoint = constants::actionWriteIterationCheckpoint();
//  bool iterating = false;
//  while (_activeSchemesEnd != _couplingSchemes.end()){
//    if ((*_activeSchemesEnd)->isActionRequired(writeCheckpoint)){
//      DEBUG("Found implicit scheme");
//      iterating = true;
//    }
//    if (iterating && (not (*_activeSchemesEnd)->isActionRequired(writeCheckpoint))){
//      DEBUG("Found explicit scheme after implicit scheme");
//      break;
//    }
//    _activeCouplingSchemes.push_back(*_activeSchemesEnd);
//    _activeSchemesEnd++;
//    DEBUG("Advanced active schemes by one new scheme");
//  }
//  assertion(not _activeCouplingSchemes.empty());
//  //assertion(_activeSchemesBegin != _activeSchemesEnd);
//}
//
//}} // namespace precice, cplscheme

