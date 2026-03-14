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
/**
 * CHECK TIME WINDOW COMPATIBILITY
 *
 * Verifies that explicit and implicit coupling schemes use compatible time window sizes.
 * In a compositional setup, the implicit scheme (e.g., FSI with Aitken acceleration) and
 * explicit schemes (e.g., sequential coupling) must coordinate their time steps.
 *
 * Requirement: explicit time window size must be an integer multiple of implicit time window size.
 * Example: implicit Dt=0.01 with explicit Dt=0.03 means explicit advances 3x per implicit window.
 *
 * This ensures synchronization points where data can be properly exchanged between schemes.
 */
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
  /**
   * REGISTER A COUPLING SUB-SCHEME
   *
   * The compositional scheme maintains exactly ONE implicit sub-scheme and any number
   * of explicit sub-schemes.  Constraints:
   *   - Only one implicit scheme is allowed (PRECICE_ASSERT guards this).
   *   - Explicit and implicit time window sizes must be compatible immediately
   *     upon registration so that mismatches are caught early.
   *
   * Registration order doesn't matter: whether explicit schemes are registered
   * before or after the implicit scheme, compatibility is checked both ways.
   */
  PRECICE_TRACE();
  if (!couplingScheme->isImplicitCouplingScheme()) {
    // Register explicit scheme and verify its time window aligns with the
    // implicit scheme (if one is already registered).
    _explicitSchemes.emplace_back(couplingScheme);

    if (_implicitScheme) {
      checkCompatibleTimeWindowSizes(*_implicitScheme, *couplingScheme);
    }
  } else {
    // Register the implicit scheme. Ensure no duplicate implicit scheme.
    PRECICE_ASSERT(_implicitScheme == nullptr);
    _implicitScheme = couplingScheme;

    // Verify compatibility against all previously registered explicit schemes.
    for (const auto &scheme : _explicitSchemes) {
      checkCompatibleTimeWindowSizes(*_implicitScheme, *scheme);
    }
  }
}

void CompositionalCouplingScheme::initialize()
{
  /**
   * INITIALIZE ALL SUB-SCHEMES AND SET UP THE ACTIVE SCHEME LIST
   *
   * On startup, ALL sub-schemes (implicit and explicit) are both initialized
   * and placed into _activeSchemes. They advance together during the first
   * time window iteration. The implicit scheme may later put explicit schemes
   * on hold (via _explicitOnHold) if it needs to iterate for convergence.
   */
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
  /**
   * Returns true if ANY sub-scheme sends initialized data.
   * A compositional scheme is considered to send initialized data when at least
   * one of its constituent schemes does (e.g., the implicit scheme initializes
   * data while explicit schemes do not).
   */
  PRECICE_TRACE();
  auto schemes              = allSchemes();
  bool sendsInitializedData = std::any_of(schemes.begin(), schemes.end(), std::mem_fn(&CouplingScheme::sendsInitializedData));
  PRECICE_DEBUG("return {}", sendsInitializedData);
  return sendsInitializedData;
}

bool CompositionalCouplingScheme::isInitialized() const
{
  /**
   * Returns true if ANY sub-scheme has been initialized.
   * Uses short-circuit OR logic: the first initialized scheme is sufficient
   * to consider the composite scheme initialized.
   */
  PRECICE_TRACE();
  auto schemes       = allSchemes();
  bool isInitialized = std::any_of(schemes.begin(), schemes.end(), std::mem_fn(&CouplingScheme::isInitialized));
  return isInitialized;
}

bool CompositionalCouplingScheme::addComputedTime(double timeToAdd)
{
  /**
   * ORCHESTRATED TIME STEPPING FOR MIXED COUPLING SCHEMES
   *
   * This is the core algorithm for managing composite coupling with implicit and explicit schemes.
   * It orchestrates how time steps are divided between:
   * - IMPLICIT scheme: iterative scheme that converges within a time window (e.g., FSI with Aitken)
   * - EXPLICIT schemes: sequential/parallel schemes that execute one time step (e.g., serial coupling)
   *
   * Key concept: SYNCHRONIZATION POINTS
   * When implicit scheme reaches convergence (end of time window), explicit schemes must pause
   * and wait. Once implicit converges, explicit schemes "catch up" to the same time.
   *
   * State machine:
   * 1. EXPLICIT PHASE (first iteration): Both implicit and explicit advance together
   * 2. IMPLICIT ITERATION: When implicit reaches time window end, it repeats. Explicit holds.
   * 3. CONVERGENCE: Implicit converges (times match). Explicit catches up.
   * 4. SYNC: All schemes synchronized at same time.
   *
   * Return: true if ANY scheme reaches end of time window in this step.
   */
  PRECICE_TRACE(timeToAdd);

  // Case 1: NO IMPLICIT SCHEME (pure explicit coupling)
  // Just advance all explicit schemes together and report when any reaches window end
  if (!isImplicitCouplingScheme()) {
    auto explicitAtWindowEnd = false;
    for (auto &scheme : _explicitSchemes) {
      explicitAtWindowEnd |= scheme->addComputedTime(timeToAdd);
    }
    return explicitAtWindowEnd;
  }

  // Case 2: WITH IMPLICIT SCHEME (composite coupling)
  // Step 1: Advance the implicit scheme (may iterate multiple times in this window)
  auto implicitAtWindowEnd = _implicitScheme->addComputedTime(timeToAdd);

  // Step 2: Check if explicit schemes are currently on hold (implicit is iterating)
  // If so, we don't advance them. Implicit iteration is more important.
  if (_explicitOnHold) {
    // Explicit schemes are frozen. Implicit scheme continues iterating.
    // No explicit schemes run in later iterations
    PRECICE_DEBUG("Explicit schemes are still on hold");
    return implicitAtWindowEnd;  // Report only implicit scheme's status
  }

  // Step 3: FIRST ITERATION - both schemes advance in parallel
  // In the first iteration through a time window, implicit and explicit both advance.
  // This provides the initial guess for implicit convergence.
  auto explicitAtWindowEnd = false;
  for (auto &scheme : _explicitSchemes) {
    explicitAtWindowEnd |= scheme->addComputedTime(timeToAdd);
  }

  // Step 4: CHECK FOR IMPLICIT CONVERGENCE
  // If implicit scheme just reached end of its time window, put explicit schemes on hold.
  // The implicit scheme will now iterate to convergence.
  if (implicitAtWindowEnd) {
    PRECICE_DEBUG("Implicit scheme reached the end of the first iteration at t={}. "
                  "Explicit schemes are on hold until convergence achieved.",
                  _implicitScheme->getTime());
    // Explicit schemes are now on hold until the implicit scheme converges
    _explicitOnHold = true;
    updateActiveSchemes();  // Update which schemes are actively participating
  }

  // Report if ANY scheme reached window end (either implicit or explicit)
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
  /**
   * FIRST DATA EXCHANGE FOR ALL ACTIVE SUB-SCHEMES
   *
   * The exchange order follows the implicit-first rule:
   *   1. If there is an implicit scheme, it always exchanges first, regardless
   *      of whether explicit schemes are on hold.
   *   2. If explicit schemes are currently on hold (implicit is still iterating),
   *      we return immediately after the implicit exchange — explicit schemes
   *      must not exchange until convergence is achieved.
   *   3. Otherwise (first iteration through the window), both implicit and
   *      explicit schemes exchange in sequence.
   */
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
  /**
   * COMPOSITE COUPLING SYNCHRONIZATION POINT
   *
   * This is the second data exchange phase where convergence is checked and schemes are
   * synchronized. The logic handles three distinct phases of a composite coupling iteration:
   *
   * PHASE 1: IMPLICIT COMPLETION
   * - Implicit scheme completes its time window (may involve multiple sub-iterations)
   * - Check if convergence has been achieved
   *
   * PHASE 2: CONVERGENCE CHECK
   * - Compare implicit and explicit times:
   *   - If implicit_time < explicit_time: NOT CONVERGED, implicit iterates, explicit waits
   *   - If implicit_time == explicit_time: CONVERGED, explicit can catch up
   *
   * PHASE 3: EXPLICIT SYNCHRONIZATION
   * - Once implicit converges, explicit schemes "fast-forward" to same time
   * - This ensures all participants synchronized before advancing to next window
   */
  PRECICE_TRACE();

  // Case 1: NO IMPLICIT SCHEME (pure explicit coupling)
  // All explicit schemes already advanced together, so just complete their exchange
  if (!isImplicitCouplingScheme()) {
    for (auto &scheme : _explicitSchemes) {
      scheme->secondExchange();
    }
    return;
  }

  /**
   * Case 2: WITH IMPLICIT SCHEME (composite coupling)
   *
   * CONVERGENCE DETECTION LOGIC:
   * Implicit scheme can lag behind explicit during iterations because when it reaches
   * time window end, it stops and starts iterating while explicit continues.
   *
   * Metaphor: Think of implicit as a runner doing laps (iterations within a window)
   * and explicit as a runner going straight. When implicit finishes a lap,
   * explicit might be ahead. Implicit then does laps until catching up.
   */

  // Step 1: Complete implicit scheme's current time window
  // This may involve sub-iterations if using iterative acceleration (e.g., Aitken)
  _implicitScheme->secondExchange();

  // Step 2: CHECK CONVERGENCE by comparing simulation times
  double implicitTime = _implicitScheme->getTime();
  double explicitTime = _explicitSchemes.front()->getTime();  // All explicit schemes in sync
  bool   iterating    = implicitTime < explicitTime;  // Implicit hasn't caught up to explicit

  if (iterating) {
    // Implicit scheme hasn't converged, implicit time is still behind explicit time
    // The implicit scheme will perform sub-iterations to bring itself forward
    // Explicit schemes must remain on hold
    PRECICE_DEBUG("Implicit scheme hasn't converged. Explicit schemes remain on hold.");
    PRECICE_ASSERT(_explicitOnHold, "Iterative scheme hasn't converged yet");
    return;  // Exit early, no need to sync explicit schemes yet
  }

  /**
   * CONVERGENCE ACHIEVED
   * Implicit scheme has caught up (implicit_time >= explicit_time).
   * Now synchronize explicit schemes to match this converged state.
   */
  PRECICE_DEBUG("Implicit scheme converged. Running explicit schemes.");

  // Step 3: SYNCHRONIZE EXPLICIT SCHEMES
  // Execute the full 3-phase cycle for explicit schemes to catch up to implicit
  // Phase 1: First exchange (receive remote data)
  for (auto &scheme : _explicitSchemes) {
    scheme->firstExchange();
  }
  // Phase 2: Second synchronization point
  for (auto &scheme : _explicitSchemes) {
    scheme->secondSynchronization();
  }
  // Phase 3: Second exchange (send computed data)
  for (auto &scheme : _explicitSchemes) {
    scheme->secondExchange();
  }

  PRECICE_DEBUG("Explicit schemes caught up. All schemes are in sync.");

  // Step 4: PREPARE FOR NEXT WINDOW
  // Explicit schemes are no longer on hold
  // Both implicit and explicit will advance together into the next time window
  _explicitOnHold = false;
  updateActiveSchemes();  // Determine which schemes participate next
}
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
  /**
   * Returns all coupling partners across every sub-scheme.
   *
   * Each sub-scheme may couple the local participant to one or more remote
   * participants. This function collects them all into a flat list, giving
   * callers a unified view of every participant this node is coupled to.
   */
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
  /**
   * Returns the local participant name.
   *
   * We use the first explicit scheme as the reference because every participant
   * in the composition shares the same local identity. If no explicit scheme
   * is registered, this case should not arise (PRECICE_ASSERT guards it).
   * The implicit scheme could also be used, but explicit is always present.
   */
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
  /**
   * COMPOSITE TIME WINDOW MANAGEMENT - SYNCHRONIZATION SCHEDULING
   *
   * In a compositional scheme with multiple coupling pairs having different time window sizes,
   * we must synchronize at all unique time windows across all couples.
   *
   * CONCRETE EXAMPLE:
   * Three participants A, B, C
   * - Couple A-B (implicit): Dt_AB = 2.0 → sync points: [2, 4, 6, 8, ...]
   * - Couple B-C (explicit): Dt_BC = 3.0 → sync points: [3, 6, 9, ...]
   *
   * Required sync points (union): [2, 3, 4, 6, 8, 9, ...]
   *
   * Time step sizes to reach consecutive sync points:
   * Starting at t=0:
   *   Step 1: advance 2.0 → reach t=2 (first sync point for A-B)
   *   Step 2: advance 1.0 → reach t=3 (first sync point for B-C)
   *   Step 3: advance 1.0 → reach t=4 (second sync point for A-B)
   *   Step 4: advance 2.0 → reach t=6 (both schemes synchronize)
   *
   * This algorithm returns the maximum step size that can be taken before hitting
   * ANY scheme's synchronization point, ensuring no missed synchronizations.
   *
   * SPECIAL CASE: When implicit scheme is iterating
   * - Explicit schemes are on hold, so only consider implicit scheme
   * - Implicit can advance its full maximum step size
   * - This avoids artificial fragmentation during implicit iterations
   */

  PRECICE_ASSERT(!_activeSchemes.empty(), "Call initialize first");

  // Return the minimum of all active schemes' max time step sizes
  // This ensures we don't overshoot any scheme's synchronization point
  double maxLength = std::transform_reduce(
      _activeSchemes.begin(), _activeSchemes.end(),
      std::numeric_limits<double>::max(),
      [](double a, double b) { return std::min(a, b); },  // Take minimum across all schemes
      std::mem_fn(&CouplingScheme::getNextTimeStepMaxSize));
  PRECICE_DEBUG("return {}", maxLength);
  return maxLength;
}

bool CompositionalCouplingScheme::isCouplingOngoing() const
{
  /**
   * Reports whether the overall coupling simulation is still running.
   *
   * When an implicit scheme is present it is always the authority: it may be
   * mid-iteration (behind the explicit time) or fully in sync, but in either
   * case its "ongoing" flag correctly reflects the global simulation state.
   * All explicit schemes are always time-synchronized once the implicit scheme
   * converges, so querying just the implicit scheme (or the first explicit
   * scheme when there is no implicit) is sufficient.
   */
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
  /**
   * MAINTAIN THE ACTIVE SCHEME SET
   *
   * _activeSchemes controls which sub-schemes participate in the current
   * time step. It is updated whenever implicit-convergence status changes:
   *
   *   Case 1: No implicit scheme
   *     → Explicit schemes are always in sync and always active.
   *       Nothing to do.
   *
   *   Case 2: Implicit scheme present, _explicitOnHold == true
   *     → Implicit scheme is iterating. Only it remains active.
   *       Explicit schemes are excluded until convergence.
   *
   *   Case 3: Implicit scheme present, _explicitOnHold == false
   *     → All schemes (implicit + explicit) are active and in sync.
   */
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
  /**
   * Returns a flat list of all registered sub-schemes (implicit + explicit).
   *
   * The implicit scheme, if present, is appended at the back of the explicit
   * list. Callers that need to visit every sub-scheme regardless of hold-state
   * (e.g., finalize(), getCouplingPartners()) use this function.
   *
   * Contrast with activeOrAllSchemes(): that function returns only the schemes
   * that are currently advancing (may exclude explicit schemes when on hold).
   */
  std::vector<CouplingScheme *> cpls(_explicitSchemes.size());
  if (_implicitScheme) {
    cpls.push_back(_implicitScheme.get());
  }
  std::transform(_explicitSchemes.begin(), _explicitSchemes.end(), cpls.begin(), std::mem_fn(&PtrCouplingScheme::get));
  return cpls;
}

std::vector<CouplingScheme *> CompositionalCouplingScheme::activeOrAllSchemes() const
{
  /**
   * Returns the currently active sub-schemes, or all schemes if none are active yet.
   *
   * _activeSchemes is managed by updateActiveSchemes() and reflects which
   * sub-schemes are advancing in the current time step:
   *   - Before initialize():  empty → falls back to allSchemes().
   *   - During implicit hold: only the implicit scheme is active.
   *   - Normal advance:       all schemes (implicit + explicit) are active.
   *
   * This function is used by methods that should only act on schemes that are
   * currently participating (e.g., isActionRequired, markActionFulfilled).
   */
  return _activeSchemes.empty() ? allSchemes() : _activeSchemes;
}

bool CompositionalCouplingScheme::isImplicitCouplingScheme() const
{
  return _implicitScheme != nullptr;
}

bool CompositionalCouplingScheme::hasConverged() const
{
  /**
   * Convergence status of the composite scheme.
   *
   * A compositional scheme without an implicit sub-scheme is always considered
   * converged (explicit coupling never iterates). When an implicit scheme is
   * present, convergence is entirely determined by whether it has converged.
   * Explicit schemes are temporarily frozen during implicit iterations, so
   * they contribute no independent convergence criterion.
   */
  if (!isImplicitCouplingScheme()) {
    return true;
  }

  return _implicitScheme->hasConverged();
}

bool CompositionalCouplingScheme::requiresSubsteps() const
{
  /**
   * Returns true if ANY sub-scheme requires substeps within a time window.
   *
   * Substeps are needed when a solver uses sub-cycling (smaller internal
   * time steps than the coupling window). If at least one constituent scheme
   * requires substeps, the composite scheme must also support them so that
   * intermediate data exchange happens at each sub-step boundary.
   */
  for (auto scheme : allSchemes()) {
    if (scheme->requiresSubsteps()) {
      return true;
    }
  }
  return false;
}

ImplicitData CompositionalCouplingScheme::implicitDataToReceive() const
{
  /**
   * Returns the implicit data that the local participant expects to receive.
   *
   * Only the implicit sub-scheme participates in iterative data exchange,
   * so only its implicit data specification is relevant here. Explicit
   * sub-schemes exchange data once per time window without iteration,
   * so they contribute nothing to implicit-receive bookkeeping.
   *
   * Returns an empty ImplicitData if no implicit scheme is registered.
   */
  if (_implicitScheme) {
    return _implicitScheme->implicitDataToReceive();
  }
  return {};
}

} // namespace precice::cplscheme
