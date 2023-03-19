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

namespace {
/// Required by STL algos, which get confused by std::min overloads
template <typename T>
T min(T a, T b)
{
  return std::min(a, b);
}
} // namespace

namespace precice::cplscheme {

void CompositionalCouplingScheme::addCouplingScheme(
    const PtrCouplingScheme &couplingScheme)
{
  PRECICE_TRACE();

  if (!couplingScheme->isImplicitCouplingScheme()) {
    _explicitSchemes.emplace_back(couplingScheme);
    return;
  }

  PRECICE_ASSERT(_implicitScheme == nullptr);
  _implicitScheme = couplingScheme;
}

void CompositionalCouplingScheme::initialize(
    double startTime,
    int    startTimeWindow)
{
  PRECICE_TRACE(startTime, startTimeWindow);
  for (const auto scheme : allSchemes()) {
    scheme->initialize(startTime, startTimeWindow);
  }
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

void CompositionalCouplingScheme::addComputedTime(double timeToAdd)
{
  PRECICE_TRACE(timeToAdd);

  for (const auto scheme : schemesToRun()) {
    scheme->addComputedTime(timeToAdd);
  }
}

CouplingScheme::ChangedMeshes CompositionalCouplingScheme::firstSynchronization(const CouplingScheme::ChangedMeshes &changes)
{
  PRECICE_TRACE();
  PRECICE_ASSERT(changes.empty());
  CouplingScheme::ChangedMeshes totalChanges;
  for (const auto scheme : schemesToRun()) {
    auto remoteChanges = scheme->firstSynchronization(changes);
    totalChanges.insert(totalChanges.end(), remoteChanges.begin(), remoteChanges.end());
  }
  return totalChanges;
}

void CompositionalCouplingScheme::firstExchange()
{
  PRECICE_TRACE();
  for (const auto scheme : schemesToRun()) {
    scheme->firstExchange();
  }
}

CouplingScheme::ChangedMeshes CompositionalCouplingScheme::secondSynchronization()
{
  PRECICE_TRACE();
  CouplingScheme::ChangedMeshes totalChanges;
  for (const auto scheme : schemesToRun()) {
    auto remoteChanges = scheme->secondSynchronization();
    totalChanges.insert(totalChanges.end(), remoteChanges.begin(), remoteChanges.end());
  }
  return totalChanges;
}

void CompositionalCouplingScheme::secondExchange()
{
  PRECICE_TRACE();
  for (const auto scheme : schemesToRun()) {
    scheme->secondExchange();
  }
  if (_implicitScheme) {
    _iterating = !_implicitScheme->hasConverged();
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
  PRECICE_TRACE();
  std::vector<std::string> partners;
  for (const auto scheme : allSchemes()) {
    auto subpartners = scheme->getCouplingPartners();
    partners.insert(partners.end(), subpartners.begin(), subpartners.end());
  }
  return partners;
}

bool CompositionalCouplingScheme::willDataBeExchanged(double lastSolverTimestepLength) const
{
  PRECICE_TRACE(lastSolverTimestepLength);
  auto schemes         = allSchemes();
  bool willBeExchanged = std::any_of(schemes.begin(), schemes.end(),
                                     [lastSolverTimestepLength](const auto &cpl) { return cpl->willDataBeExchanged(lastSolverTimestepLength); });
  PRECICE_DEBUG("return {}", willBeExchanged);
  return willBeExchanged;
}

bool CompositionalCouplingScheme::hasDataBeenReceived() const
{
  PRECICE_TRACE();

  auto schemes         = allSchemes();
  bool hasBeenReceived = std::any_of(schemes.begin(), schemes.end(), std::mem_fn(&CouplingScheme::hasDataBeenReceived));
  PRECICE_DEBUG("return {}", hasBeenReceived);
  return hasBeenReceived;
}

void CompositionalCouplingScheme::storeReceiveData(double relativeDt)
{
  PRECICE_TRACE();
  for (auto scheme : allSchemes()) {
    scheme->storeReceiveData(relativeDt);
  }
}

void CompositionalCouplingScheme::loadReceiveDataFromStorage(double relativeDt)
{
  PRECICE_TRACE();
  for (auto scheme : allSchemes()) {
    scheme->loadReceiveDataFromStorage(relativeDt);
  }
}

double CompositionalCouplingScheme::getTime() const
{
  PRECICE_TRACE();
  auto schemes = allSchemes();
  auto time    = std::transform_reduce(
      schemes.begin(), schemes.end(),
      std::numeric_limits<double>::max(),
      ::min<double>,
      std::mem_fn(&CouplingScheme::getTime));
  PRECICE_DEBUG("return {}", time);
  return time;
}

int CompositionalCouplingScheme::getTimeWindows() const
{
  PRECICE_TRACE();
  auto schemes     = allSchemes();
  auto timeWindows = std::transform_reduce(
      schemes.begin(), schemes.end(),
      std::numeric_limits<int>::max(),
      ::min<int>,
      std::mem_fn(&CouplingScheme::getTimeWindows));
  PRECICE_DEBUG("return {}", timeWindows);
  return timeWindows;
}

bool CompositionalCouplingScheme::hasTimeWindowSize() const
{
  PRECICE_TRACE();
  auto schemes = allSchemes();
  bool hasIt   = std::any_of(schemes.begin(), schemes.end(), std::mem_fn(&CouplingScheme::hasTimeWindowSize));
  PRECICE_DEBUG("return {}", hasIt);
  return hasIt;
}

double CompositionalCouplingScheme::getTimeWindowSize() const
{
  PRECICE_TRACE();
  auto   schemes        = allSchemes();
  double timeWindowSize = std::transform_reduce(
      schemes.begin(), schemes.end(), std::numeric_limits<double>::max(),
      ::min<double>,
      std::mem_fn(&CouplingScheme::getTimeWindowSize));
  PRECICE_DEBUG("return {}", timeWindowSize);
  return timeWindowSize;
}

double CompositionalCouplingScheme::getNextTimestepMaxLength() const
{
  PRECICE_TRACE();
  auto   schemes   = allSchemes();
  double maxLength = std::transform_reduce(
      schemes.begin(), schemes.end(), std::numeric_limits<double>::max(),
      ::min<double>,
      std::mem_fn(&CouplingScheme::getNextTimestepMaxLength));
  PRECICE_DEBUG("return {}", maxLength);
  return maxLength;
}

bool CompositionalCouplingScheme::isCouplingOngoing() const
{
  PRECICE_TRACE();
  auto schemes   = allSchemes();
  bool isOngoing = std::any_of(schemes.begin(), schemes.end(), std::mem_fn(&CouplingScheme::isCouplingOngoing));
  PRECICE_DEBUG("return {}", isOngoing);
  return isOngoing;
}

bool CompositionalCouplingScheme::isTimeWindowComplete() const
{
  PRECICE_TRACE();
  auto schemes    = allSchemes();
  bool isComplete = std::all_of(schemes.begin(), schemes.end(), std::mem_fn(&CouplingScheme::isTimeWindowComplete));
  PRECICE_DEBUG("return {}", isComplete);
  return isComplete;
}

bool CompositionalCouplingScheme::isActionRequired(
    Action action) const
{
  PRECICE_TRACE();
  bool isRequired = false;
  for (auto scheme : allSchemes()) {
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
  for (auto scheme : allSchemes()) {
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
  for (auto scheme : allSchemes()) {
    if (scheme->isActionRequired(action)) {
      scheme->markActionFulfilled(action);
    }
  }
}

void CompositionalCouplingScheme::requireAction(
    Action action)
{
  PRECICE_TRACE();
  for (auto scheme : schemesToRun()) {
    scheme->requireAction(action);
  }
}

std::string CompositionalCouplingScheme::printCouplingState() const
{
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

std::vector<CouplingScheme *> CompositionalCouplingScheme::schemesToRun() const
{
  if (_iterating) {
    PRECICE_DEBUG("Rerunning implicit scheme");
    return {_implicitScheme.get()};
  }
  PRECICE_DEBUG("Running all schemes");
  return allSchemes();
}

std::vector<CouplingScheme *> CompositionalCouplingScheme::allSchemes() const
{
  std::vector<CouplingScheme *> cpls(_explicitSchemes.size());
  std::transform(_explicitSchemes.begin(), _explicitSchemes.end(), cpls.begin(), std::mem_fn(&PtrCouplingScheme::get));
  if (_implicitScheme) {
    cpls.push_back(_implicitScheme.get());
  }
  return cpls;
}

bool CompositionalCouplingScheme::isImplicitCouplingScheme() const
{
  return _implicitScheme != nullptr;
}

bool CompositionalCouplingScheme::hasConverged() const
{
  if (!_implicitScheme) {
    return true;
  }

  return _implicitScheme->hasConverged();
}

} // namespace precice::cplscheme
