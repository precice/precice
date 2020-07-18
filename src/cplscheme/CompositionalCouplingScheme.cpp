#include "CompositionalCouplingScheme.hpp"
#include <algorithm>
#include <limits>
#include <memory>
#include <ostream>
#include "Constants.hpp"
#include "cplscheme/CouplingScheme.hpp"
#include "cplscheme/SharedPointer.hpp"
#include "logging/LogMacros.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace cplscheme {

void CompositionalCouplingScheme::addCouplingScheme(
    PtrCouplingScheme couplingScheme)
{
  PRECICE_TRACE();
  _couplingSchemes.emplace_back(couplingScheme);
}

void CompositionalCouplingScheme::initialize(
    double startTime,
    int    startTimeWindow)
{
  PRECICE_TRACE(startTime, startTimeWindow);
  for (Scheme scheme : _couplingSchemes) {
    scheme.scheme->initialize(startTime, startTimeWindow);
  }
  determineActiveCouplingSchemes();
}

bool CompositionalCouplingScheme::isInitialized() const
{
  PRECICE_TRACE();
  bool isInitialized = true;
  for (Scheme scheme : _couplingSchemes) {
    isInitialized &= scheme.scheme->isInitialized();
  }
  PRECICE_DEBUG("return " << isInitialized);
  return isInitialized;
}

bool CompositionalCouplingScheme::sendsInitializedData() const
{
  PRECICE_TRACE();
  bool sendsInitializedData = false;
  for (Scheme scheme : _couplingSchemes) {
    sendsInitializedData |= scheme.scheme->sendsInitializedData();
  }
  PRECICE_DEBUG("return " << sendsInitializedData);
  return sendsInitializedData;
}

bool CompositionalCouplingScheme::receivesInitializedData() const
{
  PRECICE_TRACE();
  bool receivesInitializedData = false;
  for (Scheme scheme : _couplingSchemes) {
    receivesInitializedData |= scheme.scheme->receivesInitializedData();
  }
  PRECICE_DEBUG("return " << receivesInitializedData);
  return receivesInitializedData;
}

void CompositionalCouplingScheme::initializeData()
{
  PRECICE_TRACE();
  for (Scheme scheme : _couplingSchemes) {
    scheme.scheme->initializeData();
  }
}

void CompositionalCouplingScheme::addComputedTime(double timeToAdd)
{
  PRECICE_TRACE(timeToAdd);
  for (SchemesIt it = _activeSchemesBegin; it != _activeSchemesEnd; it++) {
    if (not it->onHold) {
      it->scheme->addComputedTime(timeToAdd);
    }
  }
  _lastAddedTime += timeToAdd;
}

void CompositionalCouplingScheme::advance()
{
  PRECICE_TRACE();
  bool moreSchemesToHandle = false;
  do {
    for (SchemesIt it = _activeSchemesBegin; it != _activeSchemesEnd; it++) {
      if (not it->onHold) {
        it->scheme->advance();
      }
    }
    moreSchemesToHandle = determineActiveCouplingSchemes();
    if (moreSchemesToHandle) {
      // The new schemes to be handled in this advance also need the time that
      // has been computed so far. This time can't be added in the solver call
      // to addComputedTime(), since there the schemes are not active yet.
      addComputedTime(_lastAddedTime);
    }
  } while (moreSchemesToHandle);
  _lastAddedTime = 0.0;
}

void CompositionalCouplingScheme::finalize()
{
  PRECICE_TRACE();
  for (Scheme scheme : _couplingSchemes) {
    scheme.scheme->finalize();
  }
}

std::vector<std::string> CompositionalCouplingScheme::getCouplingPartners() const
{
  PRECICE_TRACE();
  std::vector<std::string> partners;
  std::vector<std::string> subpartners;
  for (Scheme scheme : _couplingSchemes) {
    subpartners = scheme.scheme->getCouplingPartners();
    partners.insert(partners.end(), subpartners.begin(), subpartners.end());
  }
  return partners;
}

bool CompositionalCouplingScheme::willDataBeExchanged(double lastSolverTimestepLength) const
{
  PRECICE_TRACE(lastSolverTimestepLength);
  bool willBeExchanged = false;
  for (SchemesIt it = _activeSchemesBegin; it != _activeSchemesEnd; it++) {
    if (not it->onHold) {
      willBeExchanged |= it->scheme->willDataBeExchanged(lastSolverTimestepLength);
    }
  }
  PRECICE_DEBUG("return " << willBeExchanged);
  return willBeExchanged;
}

bool CompositionalCouplingScheme::hasDataBeenReceived() const
{
  PRECICE_TRACE();
  bool hasBeenReceived = false;
  // Question: Does it suffice to only check the active ones?
  for (SchemesIt it = _activeSchemesBegin; it != _activeSchemesEnd; it++) {
    if (not it->onHold) {
      hasBeenReceived |= it->scheme->hasDataBeenReceived();
    }
  }
  PRECICE_DEBUG("return " << hasBeenReceived);
  return hasBeenReceived;
}

double CompositionalCouplingScheme::getTime() const
{
  PRECICE_TRACE();
  double time = std::numeric_limits<double>::max();
  for (Scheme scheme : _couplingSchemes) {
    if (not scheme.onHold) {
      time = std::min(time, scheme.scheme->getTime());
    }
  }
  PRECICE_DEBUG("return " << time);
  return time;
}

int CompositionalCouplingScheme::getTimeWindows() const
{
  PRECICE_TRACE();
  int timeWindows = std::numeric_limits<int>::max();
  for (Scheme scheme : _couplingSchemes) {
    if (not scheme.onHold) {
      timeWindows = std::min(timeWindows, scheme.scheme->getTimeWindows());
    }
  }
  PRECICE_DEBUG("return " << timeWindows);
  return timeWindows;
}

bool CompositionalCouplingScheme::hasTimeWindowSize() const
{
  PRECICE_TRACE();
  bool hasIt = false;
  for (Scheme scheme : _couplingSchemes) {
    hasIt |= scheme.scheme->hasTimeWindowSize();
  }
  PRECICE_DEBUG("return " << hasIt);
  return hasIt;
}

double CompositionalCouplingScheme::getTimeWindowSize() const
{
  PRECICE_TRACE();
  double timeWindowSize = std::numeric_limits<double>::max();
  for (Scheme scheme : _couplingSchemes) {
    if (scheme.scheme->getTimeWindowSize() < timeWindowSize) {
      timeWindowSize = scheme.scheme->getTimeWindowSize();
    }
  }
  PRECICE_DEBUG("return " << timeWindowSize);
  return timeWindowSize;
}

double CompositionalCouplingScheme::getThisTimeWindowRemainder() const
{
  PRECICE_TRACE();
  double maxRemainder = 0.0;
  for (Scheme scheme : _couplingSchemes) {
    if (not scheme.onHold) {
      if (scheme.scheme->getThisTimeWindowRemainder() > maxRemainder) {
        maxRemainder = scheme.scheme->getThisTimeWindowRemainder();
      }
    }
  }
  PRECICE_DEBUG("return " << maxRemainder);
  return maxRemainder;
}

double CompositionalCouplingScheme::getNextTimestepMaxLength() const
{
  PRECICE_TRACE();
  double maxLength = std::numeric_limits<double>::max();
  for (Scheme scheme : _couplingSchemes) {
    if (not scheme.onHold) {
      maxLength = std::min(maxLength, scheme.scheme->getNextTimestepMaxLength());
    }
  }
  PRECICE_DEBUG("return " << maxLength);
  return maxLength;
}

bool CompositionalCouplingScheme::isCouplingOngoing() const
{
  PRECICE_TRACE();
  bool isOngoing = false;
  for (Scheme scheme : _couplingSchemes) {
    isOngoing |= scheme.scheme->isCouplingOngoing();
  }
  PRECICE_DEBUG("return " << isOngoing);
  return isOngoing;
}

bool CompositionalCouplingScheme::isTimeWindowComplete() const
{
  PRECICE_TRACE();
  bool isComplete = true;
  for (Scheme scheme : _couplingSchemes) {
    isComplete &= scheme.scheme->isTimeWindowComplete();
  }
  PRECICE_DEBUG("return " << isComplete);
  return isComplete;
}

bool CompositionalCouplingScheme::isActionRequired(
    const std::string &actionName) const
{
  PRECICE_TRACE(actionName);
  bool isRequired = false;
  for (Scheme scheme : _couplingSchemes) {
    if (not scheme.onHold) {
      isRequired |= scheme.scheme->isActionRequired(actionName);
    }
  }
  PRECICE_DEBUG("return " << isRequired);
  return isRequired;
}

void CompositionalCouplingScheme::markActionFulfilled(
    const std::string &actionName)
{
  PRECICE_TRACE(actionName);
  for (Scheme scheme : _couplingSchemes) {
    if (not scheme.onHold) {
      scheme.scheme->markActionFulfilled(actionName);
    }
  }
}

void CompositionalCouplingScheme::requireAction(
    const std::string &actionName)
{
  PRECICE_TRACE(actionName);
  for (Scheme scheme : _couplingSchemes) {
    scheme.scheme->requireAction(actionName);
  }
}

std::string CompositionalCouplingScheme::printCouplingState() const
{
  std::string              state;
  std::vector<std::string> partners;
  for (Scheme scheme : _couplingSchemes) {
    if (not state.empty()) {
      state += "\n";
    }
    partners = scheme.scheme->getCouplingPartners();
    state += partners[0];
    state += ": ";
    state += scheme.scheme->printCouplingState();
  }
  return state;
}

bool CompositionalCouplingScheme::determineActiveCouplingSchemes()
{
  PRECICE_TRACE();
  bool        newActiveSchemes = false;
  std::string writeCheckpoint  = constants::actionWriteIterationCheckpoint();
  std::string readCheckpoint   = constants::actionReadIterationCheckpoint();
  if (_activeSchemesBegin == _activeSchemesEnd) {
    PRECICE_DEBUG("Case After Init");
    // First call after initialization of all coupling schemes. All coupling
    // schemes are set active up to (but not including) the first explicit
    // scheme after an implicit scheme.
    _activeSchemesBegin = _couplingSchemes.begin();
    _activeSchemesEnd   = _couplingSchemes.begin();
    advanceActiveCouplingSchemes();
    newActiveSchemes = true;
  } else {
    PRECICE_DEBUG("Normal Case");
    // Redetermine active schemes. First, all preceding explicit schemes are
    // removed. Then, all remaining implicit schemes are checked for the
    // convergence of iterations (this is given when an iteration checkpoint
    // should be created). If all are converged, a next set of active schemes
    // is determined.

    // Remove preceding explicit schemes
    while (_activeSchemesBegin != _activeSchemesEnd) {
      bool explicitScheme = true;
      explicitScheme &= not _activeSchemesBegin->scheme->isActionRequired(writeCheckpoint);
      explicitScheme &= not _activeSchemesBegin->scheme->isActionRequired(readCheckpoint);
      if (explicitScheme) {
        _activeSchemesBegin++;
        PRECICE_DEBUG("Remove preceding explicit scheme");
      } else {
        break;
      }
    }

    // Check implicit schemes for convergence and remove if converged
    bool converged = true;
    for (SchemesIt it = _activeSchemesBegin; it != _activeSchemesEnd; it++) {
      if (it->scheme->isActionRequired(readCheckpoint)) {
        converged = false;
        PRECICE_DEBUG("Non converged implicit scheme");
      } else if (it->scheme->isActionRequired(writeCheckpoint) || not it->scheme->isCouplingOngoing()) {
        it->onHold = true;
        PRECICE_DEBUG("Put converged/finished implicit scheme on hold");
      }
    }
    if (converged) {
      PRECICE_DEBUG("Active implicit schemes converged");
      for (SchemesIt it = _activeSchemesBegin; it != _activeSchemesEnd; it++) {
        it->onHold = false;
      }
      _activeSchemesBegin = _activeSchemesEnd;
    }

    // Determine next set of active schemes if current is empty
    if (_activeSchemesBegin == _activeSchemesEnd) {
      if (_activeSchemesBegin == _couplingSchemes.end()) {
        PRECICE_DEBUG("Through with all coupling schemes");
        // All coupling schemes are through
        _activeSchemesBegin = _couplingSchemes.begin();
        _activeSchemesEnd   = _couplingSchemes.begin();
        advanceActiveCouplingSchemes();
        // newActiveSchemes stays false, since the current it/dt is complete
      } else {
        PRECICE_DEBUG("Coupling schemes remaining");
        advanceActiveCouplingSchemes();
        newActiveSchemes = true;
      }
    }
  }
  PRECICE_DEBUG("return newActiveSchemes=" << newActiveSchemes);
  return newActiveSchemes;
}

void CompositionalCouplingScheme::advanceActiveCouplingSchemes()
{
  PRECICE_TRACE();
  std::string writeCheckpoint = constants::actionWriteIterationCheckpoint();
  bool        iterating       = false;
  while (_activeSchemesEnd != _couplingSchemes.end()) {
    if (_activeSchemesEnd->scheme->isActionRequired(writeCheckpoint)) {
      PRECICE_DEBUG("Found implicit scheme");
      iterating = true;
    }
    if (iterating && (not _activeSchemesEnd->scheme->isActionRequired(writeCheckpoint))) {
      PRECICE_DEBUG("Found explicit scheme after implicit scheme");
      break;
    }
    _activeSchemesEnd++;
    PRECICE_DEBUG("Advanced active schemes by one new scheme");
  }
  PRECICE_ASSERT(_activeSchemesBegin != _activeSchemesEnd);
}

} // namespace cplscheme
} // namespace precice
