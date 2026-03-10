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
// HELPER: Time window size compatibility check
// ==============================================
// Validates that explicit and implicit coupling schemes have compatible time windows.
// For stable convergence in multi-coupling, explicit time steps must align with implicit steps.
//
// Rule: edt (explicit) must be an integer multiple of idt (implicit)
// Example: If implicit uses dt=0.01, explicit can use 0.01, 0.02, 0.03, etc.
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
  // MULTI-COUPLING REGISTRATION
  // ============================
  // Compositional coupling combines implicit and explicit schemes:
  //   - Implicit (at most 1): Always runs, drives convergence iterations
  //   - Explicit (multiple): Second-order accurate, runs after implicit converges
  //
  // Example configuration:
  //   Participant A <-> B (implicit, convergence loop)
  //   Participant A <-> C (explicit, runs once per implicit time window)
  //
  // Validation:
  //   - Time window sizes must be compatible (explicit = n * implicit)
  //   - At most one implicit scheme allowed
  
  if (!couplingScheme->isImplicitCouplingScheme()) {
    // Register as explicit coupling (sequential, second-order)
    _explicitSchemes.emplace_back(couplingScheme);

    if (_implicitScheme) {
      checkCompatibleTimeWindowSizes(*_implicitScheme, *couplingScheme);
    }
  } else {
    // Register as implicit coupling (iterative, drives convergence)
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
  // INITIALIZATION OF MULTI-COUPLING SCHEMES
  // =========================================
  // Initialize all registered coupling schemes (implicit + explicit).
  // Both scheme types are active at startup; explicit schemes may be put
  // on hold later during implicit convergence iterations.
  
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
  // TIME WINDOW ADVANCEMENT AND IMPLICIT-EXPLICIT ORCHESTRATION
  // =============================================================
  // Manages time stepping in multi-coupled scenarios:
  //
  // CASE 1: Pure explicit coupling (no implicit scheme)
  //   - All schemes advance time in lockstep
  //   - Return true when any explicit scheme reaches window end
  //
  // CASE 2: Multi-coupling (implicit + explicit)
  //   ITERATION 1 (first):
  //     - Implicit and explicit advance together
  //     - When implicit reaches window end → explicit schemes go ON HOLD
  //
  //   ITERATIONS 2+ (implicit convergence loop):
  //     - Only implicit continues iterating (explicit frozen)
  //     - After implicit converges → explicit schemes resume

  if (!isImplicitCouplingScheme()) {
    // Pure explicit coupling: all advance together
    auto explicitAtWindowEnd = false;
    for (auto &scheme : _explicitSchemes) {
      explicitAtWindowEnd |= scheme->addComputedTime(timeToAdd);
    }
    return explicitAtWindowEnd;
  }

  // Multi-coupling case: coordinate implicit-explicit interaction
  // Step 1: Advance implicit scheme (drives iterations)
  auto implicitAtWindowEnd = _implicitScheme->addComputedTime(timeToAdd);

  // Step 2: If explicit schemes are on hold, skip them (implicit iterating)
  if (_explicitOnHold) {
    PRECICE_DEBUG("Explicit schemes are still on hold");
    return implicitAtWindowEnd;
  }

  // Step 3: In first iteration, explicit schemes advance with implicit
  auto explicitAtWindowEnd = false;
  for (auto &scheme : _explicitSchemes) {
    explicitAtWindowEnd |= scheme->addComputedTime(timeToAdd);
  }

  // Step 4: When implicit reaches window end, explicit must freeze
  if (implicitAtWindowEnd) {
    PRECICE_DEBUG("Implicit scheme reached the end of the first iteration at t={}. "
                  "Explicit schemes are on hold until convergence achieved.",
                  _implicitScheme->getTime());
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
  // FIRST DATA EXCHANGE IN MULTI-COUPLING
  // ======================================
  // Initial data exchange phase. Behavior depends on coupling configuration:
  //
  // Pure explicit: Exchange with all schemes
  // Multi-coupling: 
  //   - Iteration 1: Implicit exchanges, then explicit (if not on hold)
  //   - Iterations 2+: Only implicit exchanges (explicit frozen)

  if (!isImplicitCouplingScheme()) {
    // Pure explicit: exchange with all
    for (auto &scheme : _explicitSchemes) {
      scheme->firstExchange();
    }
    return;
  }

  // Multi-coupling: implicit always exchanges
  _implicitScheme->firstExchange();
  
  // Explicit only exchanges if not iterating implicit convergence
  if (_explicitOnHold) {
    return;
  }

  // First iteration: explicit exchanges in parallel with implicit
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
  // SECOND DATA EXCHANGE AND CONVERGENCE SYNCHRONIZATION
  // ======================================================
  // Handles convergence detection and explicit scheme activation.
  //
  // Multi-coupling flow:
  //   1. Implicit completes time window → check convergence
  //   2. If NOT converged (iterating):
  //      - Explicit remains on hold
  //      - Loop continues with implicit scheme
  //   3. If converged:
  //      - Run explicit schemes (they were frozen)
  //      - Bring explicit up to same time as implicit
  //      - Synchronize all schemes for next time window

  if (!isImplicitCouplingScheme()) {
    // Pure explicit: all schemes exchange
    for (auto &scheme : _explicitSchemes) {
      scheme->secondExchange();
    }
    return;
  }

  // Multi-coupling: Check if implicit scheme has converged
  _implicitScheme->secondExchange();

  // Convergence criterion: implicit time >= explicit time
  // (implicit has dominated the coupling loop)
  double implicitTime = _implicitScheme->getTime();
  double explicitTime = _explicitSchemes.front()->getTime();
  bool   iterating    = implicitTime < explicitTime;

  if (iterating) {
    // Implicit still converging: explicit remains frozen
    PRECICE_DEBUG("Implicit scheme hasn't converged. Explicit schemes remain on hold.");
    PRECICE_ASSERT(_explicitOnHold, "Iterative scheme hasn't converged yet");
    return;
  }
  
  // Implicit converged: release explicit schemes to catch up
  PRECICE_DEBUG("Implicit scheme converged. Running explicit schemes.");

  // Run all exchange phases for explicit schemes to complete their time window
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

  // Ready for next time window
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
  std::vector<std::string> states;
  for (const auto scheme : allSchemes()) {
    states.push_back(fmt::format("partner: {}, {}", fmt::join(scheme->getCouplingPartners(), " & "), scheme->printCouplingState()));
  }
  return fmt::format("{}", fmt::join(states, "; "));
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
  std::vector<CouplingScheme *> cpls;
  cpls.reserve(_explicitSchemes.size() + (_implicitScheme ? 1 : 0));
  
  // Add explicit schemes
  std::transform(_explicitSchemes.begin(), _explicitSchemes.end(), std::back_inserter(cpls), 
                 std::mem_fn(&PtrCouplingScheme::get));
  
  // Add implicit scheme if it exists
  if (_implicitScheme) {
    cpls.push_back(_implicitScheme.get());
  }
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
