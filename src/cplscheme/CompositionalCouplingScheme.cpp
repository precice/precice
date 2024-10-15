#include <algorithm>
#include <functional>
#include <limits>
#include <memory>
#include <numeric>
#include <ostream>

#include "CompositionalCouplingScheme.hpp"
#include "Constants.hpp"
#include "cplscheme/BaseCouplingScheme.hpp"
#include "cplscheme/CouplingScheme.hpp"
#include "cplscheme/SharedPointer.hpp"
#include "cplscheme/tests/DummyCouplingScheme.hpp"
#include "logging/LogMacros.hpp"
#include "utils/assertion.hpp"

namespace precice::cplscheme {

namespace {
bool compatibleTimeWindowSizes(const CouplingScheme &impl, const CouplingScheme &expl)
{
  if (!impl.hasTimeWindowSize() || !expl.hasTimeWindowSize()) {
    return true;
  }
  double idt = impl.getTimeWindowSize();
  double edt = expl.getTimeWindowSize();
  // edt needs to be an integer multiple of idt
  return math::equals(edt, idt) || math::equals(std::remainder(edt, idt), 0.0);
}
} // namespace

void CompositionalCouplingScheme::checkCompatibleTimeWindowSizes(const CouplingScheme &impl, const CouplingScheme &expl) const
{
  if (compatibleTimeWindowSizes(impl, expl)) {
    return;
  }
  auto   local        = expl.localParticipant();
  auto   explPartners = expl.getCouplingPartners();
  auto   implPartners = impl.getCouplingPartners();
  double edt          = expl.getTimeWindowSize();
  double idt          = impl.getTimeWindowSize();

  PRECICE_WARN(
      "The participant {} is implicitly coupled to {} with time-window-size {} and explicitly coupled to {} with time-window-size {}, which isn't well-defined. "
      "Explicit time windows should be aligned with implicit time windows. "
      "Choose time window sizes so that the explicit (currently {}) is an integer multiple of the implicit (currently {}).",
      local, implPartners.front(), idt, explPartners.front(), edt,
      edt, idt);
}

void CompositionalCouplingScheme::addCouplingScheme(
    const PtrCouplingScheme &couplingScheme)
{
  PRECICE_TRACE();
  if (!couplingScheme->isImplicitCouplingScheme()) {
    _explicitSchemes.emplace_back(couplingScheme);

    if (_implicitScheme) {
      checkCompatibleTimeWindowSizes(*_implicitScheme, *couplingScheme);
    }
  } else {
    PRECICE_ASSERT(_implicitScheme == nullptr);
    _implicitScheme = couplingScheme;

    for (const auto &scheme : _explicitSchemes) {
      checkCompatibleTimeWindowSizes(*_implicitScheme, *scheme);
    }
  }
}

void CompositionalCouplingScheme::initialize()
{
  PRECICE_TRACE();
  PRECICE_ASSERT(_activeSchemes.empty());
  for (const auto scheme : allSchemes()) {
    scheme->initialize();
    _activeSchemes.push_back(scheme);
  }
}

void CompositionalCouplingScheme::reinitialize()
{
  PRECICE_TRACE();
  PRECICE_UNREACHABLE("Not implemented and not allowed");
}

bool CompositionalCouplingScheme::sendsInitializedData() const
{
  PRECICE_TRACE();
  auto schemes              = allSchemes();
  bool sendsInitializedData = std::any_of(schemes.begin(), schemes.end(), std::mem_fn(&CouplingScheme::sendsInitializedData));
  PRECICE_DEBUG("return {}", sendsInitializedData);
  return sendsInitializedData;
}

bool CompositionalCouplingScheme::isInitialized() const
{
  PRECICE_TRACE();
  auto schemes       = allSchemes();
  bool isInitialized = std::any_of(schemes.begin(), schemes.end(), std::mem_fn(&CouplingScheme::isInitialized));
  return isInitialized;
}

bool CompositionalCouplingScheme::addComputedTime(double timeToAdd)
{
  PRECICE_TRACE(timeToAdd);

  if (!isImplicitCouplingScheme()) {
    auto explicitAtWindowEnd = false;
    for (auto &scheme : _explicitSchemes) {
      explicitAtWindowEnd |= scheme->addComputedTime(timeToAdd);
    }
    return explicitAtWindowEnd;
  }

  // Reaching the timewindow end of the implicit scheme requires explicit coupling schemes to freeze
  auto implicitAtWindowEnd = _implicitScheme->addComputedTime(timeToAdd);

  // We are done if the implicit scheme is iterating
  if (_explicitOnHold) {
    // No explicit schemes run in later iterations
    PRECICE_DEBUG("Explicit schemes are still on hold");
    return implicitAtWindowEnd;
  }

  // We are in the first iteration so explicit schemes need to step forward
  auto explicitAtWindowEnd = false;
  for (auto &scheme : _explicitSchemes) {
    explicitAtWindowEnd |= scheme->addComputedTime(timeToAdd);
  }

  if (implicitAtWindowEnd) {
    PRECICE_DEBUG("Implicit scheme reached the end of the first iteration at t={}. "
                  "Explicit schemes are on hold until convergence achieved.",
                  _implicitScheme->getTime());
    // explicit schemes are on hold until converged
    _explicitOnHold = true;
    updateActiveSchemes();
  }

  return implicitAtWindowEnd || explicitAtWindowEnd;
}

CouplingScheme::ChangedMeshes CompositionalCouplingScheme::firstSynchronization(const CouplingScheme::ChangedMeshes &changes)
{
  PRECICE_TRACE();
  PRECICE_ASSERT(changes.empty());
  PRECICE_ASSERT(!_activeSchemes.empty(), "Call initialize first");
  CouplingScheme::ChangedMeshes totalChanges;
  for (const auto scheme : _activeSchemes) {
    auto remoteChanges = scheme->firstSynchronization(changes);
    totalChanges.insert(totalChanges.end(), remoteChanges.begin(), remoteChanges.end());
  }
  return totalChanges;
}

void CompositionalCouplingScheme::firstExchange()
{
  PRECICE_TRACE();

  if (!isImplicitCouplingScheme()) {
    for (auto &scheme : _explicitSchemes) {
      scheme->firstExchange();
    }
    return;
  }

  // The implicit scheme either just reached the end of the first time window, or is already iterating.
  _implicitScheme->firstExchange();
  if (_explicitOnHold) {
    return;
  }

  // The implicit scheme hasn't reached the end of the first time window yet.
  for (auto &scheme : _explicitSchemes) {
    scheme->firstExchange();
  }
}

CouplingScheme::ChangedMeshes CompositionalCouplingScheme::secondSynchronization()
{
  PRECICE_TRACE();
  PRECICE_ASSERT(!_activeSchemes.empty(), "Call initialize first");
  CouplingScheme::ChangedMeshes totalChanges;
  for (const auto scheme : _activeSchemes) {
    auto remoteChanges = scheme->secondSynchronization();
    totalChanges.insert(totalChanges.end(), remoteChanges.begin(), remoteChanges.end());
  }
  return totalChanges;
}

void CompositionalCouplingScheme::secondExchange()
{
  PRECICE_TRACE();

  if (!isImplicitCouplingScheme()) {
    for (auto &scheme : _explicitSchemes) {
      scheme->secondExchange();
    }
    return;
  }

  // First complete the time window of the implicit scheme to determine if the scheme has converged
  _implicitScheme->secondExchange();

  double implicitTime = _implicitScheme->getTime();
  double explicitTime = _explicitSchemes.front()->getTime();
  bool   iterating    = implicitTime < explicitTime; // TODO numeric check?

  if (iterating) {
    PRECICE_DEBUG("Implicit scheme hasn't converged. Explicit schemes remain on hold.");
    PRECICE_ASSERT(_explicitOnHold, "Iterative scheme hasn't converged yet");
    return;
  }
  PRECICE_DEBUG("Implicit scheme converged. Running explicit schemes.");

  for (auto &scheme : _explicitSchemes) {
    scheme->firstExchange();
  }
  for (auto &scheme : _explicitSchemes) {
    scheme->secondSynchronization();
  }
  for (auto &scheme : _explicitSchemes) {
    scheme->secondExchange();
  }
  PRECICE_DEBUG("Explicit schemes caught up. All schemes are in sync.");

  _explicitOnHold = false;
  updateActiveSchemes();
}

void CompositionalCouplingScheme::finalize()
{
  PRECICE_TRACE();
  for (const auto scheme : allSchemes()) {
    scheme->finalize();
  }
}

std::vector<std::string> CompositionalCouplingScheme::getCouplingPartners() const
{
  PRECICE_TRACE();
  std::vector<std::string> partners;
  for (const auto scheme : allSchemes()) {
    auto subpartners = scheme->getCouplingPartners();
    partners.insert(partners.end(), subpartners.begin(), subpartners.end());
  }
  return partners;
}

std::string CompositionalCouplingScheme::localParticipant() const
{
  PRECICE_TRACE();
  PRECICE_ASSERT(!_explicitSchemes.empty());
  return _explicitSchemes.front()->localParticipant();
}

bool CompositionalCouplingScheme::willDataBeExchanged(double lastSolverTimeStepSize) const
{
  PRECICE_TRACE(lastSolverTimeStepSize);
  // @TODO
  auto schemes         = allSchemes();
  bool willBeExchanged = std::any_of(schemes.begin(), schemes.end(),
                                     [lastSolverTimeStepSize](const auto &cpl) { return cpl->willDataBeExchanged(lastSolverTimeStepSize); });
  PRECICE_DEBUG("return {}", willBeExchanged);
  return willBeExchanged;
}

bool CompositionalCouplingScheme::hasDataBeenReceived() const
{
  PRECICE_TRACE();
  // @TODO
  auto schemes         = allSchemes();
  bool hasBeenReceived = std::any_of(schemes.begin(), schemes.end(), std::mem_fn(&CouplingScheme::hasDataBeenReceived));
  PRECICE_DEBUG("return {}", hasBeenReceived);
  return hasBeenReceived;
}

double CompositionalCouplingScheme::getTime() const
{
  PRECICE_TRACE();
  if (_implicitScheme) {
    // Either the implicit scheme is behind while iterating or all schemes are in sync
    return _implicitScheme->getTime();
  } else {
    // Explicit schemes are always in sync
    return _explicitSchemes.front()->getTime();
  }
}

double CompositionalCouplingScheme::getTimeWindowStart() const
{
  PRECICE_TRACE();
  // In case of various time window sizes, use the scheme with the earliest start.
  auto schemes = allSchemes();
  return std::transform_reduce(
      schemes.begin(), schemes.end(),
      std::numeric_limits<double>::max(),
      [](double a, double b) { return std::min(a, b); },
      std::mem_fn(&CouplingScheme::getTimeWindowStart));
}

int CompositionalCouplingScheme::getTimeWindows() const
{
  PRECICE_TRACE();
  // In case of various time window sizes, use the scheme with the most timeWindows.
  auto schemes     = allSchemes();
  auto timeWindows = std::transform_reduce(
      schemes.begin(), schemes.end(),
      std::numeric_limits<int>::max(),
      [](int a, int b) { return std::min(a, b); },
      std::mem_fn(&CouplingScheme::getTimeWindows));
  PRECICE_DEBUG("return {}", timeWindows);
  return timeWindows;
}

bool CompositionalCouplingScheme::hasTimeWindowSize() const
{
  // TODO Move to BaseCouplingScheme
  PRECICE_TRACE();
  PRECICE_UNREACHABLE("This may not be called for a CompositionalCouplingScheme");
  return true;
}

double CompositionalCouplingScheme::getTimeWindowSize() const
{
  // TODO Move to BaseCouplingScheme
  PRECICE_TRACE();
  PRECICE_UNREACHABLE("This may not be called for a CompositionalCouplingScheme");
  return -1.0;
}

double CompositionalCouplingScheme::getNextTimeStepMaxSize() const
{
  PRECICE_TRACE();
  PRECICE_ASSERT(!_activeSchemes.empty(), "Call initialize first");

  /* For compositional schemes with mixed time window sizes, this is the max time step size to the next end of a time window.
   *
   * As a concrete example we have 3 solvers A,B,C which are coupled A-B with Dt=2 and B-C with Dt=3.
   *
   * Then the synchronization points for A-B are (2, 4, 6, 8, ...) and for B-C are (3, 6, 9, ...) .
   * The compositional scheme has to synchronize at the union of all these points (2, 3, 4, 6, 8, 9, ...)
   *
   * Hence the getNextTimeStepMaxSize() is the difference between these points.
   * Starting at t=0 and always advancing the full getNextTimeStepMaxSize() will lead to time steps of sizes (2,1,1,2,2,1,...)
   *
   * Exception is when the implicit scheme is iterating. In that moment the explicit schemes are on hold and the implicit scheme may advance its full max time step size.
   */

  double maxLength = std::transform_reduce(
      _activeSchemes.begin(), _activeSchemes.end(), std::numeric_limits<double>::max(),
      [](double a, double b) { return std::min(a, b); },
      std::mem_fn(&CouplingScheme::getNextTimeStepMaxSize));
  PRECICE_DEBUG("return {}", maxLength);
  return maxLength;
}

bool CompositionalCouplingScheme::isCouplingOngoing() const
{
  PRECICE_TRACE();
  if (_implicitScheme) {
    // The implicit scheme is either the one that has to catch up, or all schemes are in sync
    return _implicitScheme->isCouplingOngoing();
  }
  // All explicit schemes are always in sync
  return _explicitSchemes.front()->isCouplingOngoing();
}

bool CompositionalCouplingScheme::isTimeWindowComplete() const
{
  PRECICE_TRACE();
  PRECICE_ASSERT(!_activeSchemes.empty(), "Call initialize first");
  // Same situation as described in getNextTimeStepMaxSize(), but we are interested in reaching these time windows
  bool isOneCompleted = std::any_of(_activeSchemes.begin(), _activeSchemes.end(), std::mem_fn(&CouplingScheme::isTimeWindowComplete));
  PRECICE_DEBUG("return {}", isOneCompleted);
  return isOneCompleted;
}

bool CompositionalCouplingScheme::isActionRequired(
    Action action) const
{
  PRECICE_TRACE();
  bool isRequired = false;

  for (auto scheme : activeOrAllSchemes()) {
    if (scheme->isActionRequired(action)) {
      isRequired = true;
      break;
    }
  }
  PRECICE_DEBUG("return {}", isRequired);
  return isRequired;
}

bool CompositionalCouplingScheme::isActionFulfilled(
    Action action) const
{
  PRECICE_TRACE();
  bool isFulfilled = false;
  for (auto scheme : activeOrAllSchemes()) {
    if (scheme->isActionFulfilled(action)) {
      isFulfilled = true;
      break;
    }
  }
  PRECICE_DEBUG("return {}", isFulfilled);
  return isFulfilled;
}

void CompositionalCouplingScheme::markActionFulfilled(
    Action action)
{
  PRECICE_TRACE();
  for (auto scheme : activeOrAllSchemes()) {
    if (scheme->isActionRequired(action)) {
      scheme->markActionFulfilled(action);
    }
  }
}

void CompositionalCouplingScheme::requireAction(
    Action action)
{
  PRECICE_TRACE();
  PRECICE_UNREACHABLE("This may not be called for a CompositionalCouplingScheme");
}

std::string CompositionalCouplingScheme::printCouplingState() const
{
  // TODO is a redesign necessary?
  std::string              state;
  std::vector<std::string> partners;
  for (const auto scheme : allSchemes()) {
    if (not state.empty()) {
      state += "\n";
    }
    partners = scheme->getCouplingPartners();
    state += partners[0];
    state += ": ";
    state += scheme->printCouplingState();
  }
  return state;
}

void CompositionalCouplingScheme::updateActiveSchemes()
{
  if (_implicitScheme == nullptr) {
    // We only have explicit schemes, which are always in sync and thus always run together
    // Active schemes never change, so there is nothing to do
    PRECICE_ASSERT(_activeSchemes.size() == _explicitSchemes.size(), "Active scheme haven't been initialized correctly!");
    return;
  }

  // There is an implicit scheme, which may be iterating

  _activeSchemes.clear();

  // The implicit scheme is either in sync or has to catch up
  _activeSchemes.push_back(_implicitScheme.get());

  if (_explicitOnHold) {
    return;
  }

  // If the implicit scheme isn't currently iterating, then all schemes are in sync
  for (auto &scheme : _explicitSchemes) {
    _activeSchemes.push_back(scheme.get());
  }
}

std::vector<CouplingScheme *> CompositionalCouplingScheme::allSchemes() const
{
  std::vector<CouplingScheme *> cpls(_explicitSchemes.size());
  if (_implicitScheme) {
    cpls.push_back(_implicitScheme.get());
  }
  std::transform(_explicitSchemes.begin(), _explicitSchemes.end(), cpls.begin(), std::mem_fn(&PtrCouplingScheme::get));
  return cpls;
}

std::vector<CouplingScheme *> CompositionalCouplingScheme::activeOrAllSchemes() const
{
  return _activeSchemes.empty() ? allSchemes() : _activeSchemes;
}

bool CompositionalCouplingScheme::isImplicitCouplingScheme() const
{
  return _implicitScheme != nullptr;
}

bool CompositionalCouplingScheme::hasConverged() const
{
  if (!isImplicitCouplingScheme()) {
    return true;
  }

  return _implicitScheme->hasConverged();
}

bool CompositionalCouplingScheme::requiresSubsteps() const
{
  for (auto scheme : allSchemes()) {
    if (scheme->requiresSubsteps()) {
      return true;
    }
  }
  return false;
}

ImplicitData CompositionalCouplingScheme::implicitDataToReceive() const
{
  if (_implicitScheme) {
    return _implicitScheme->implicitDataToReceive();
  }
  return {};
}

} // namespace precice::cplscheme
