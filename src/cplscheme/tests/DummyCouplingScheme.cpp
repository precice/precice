#include "DummyCouplingScheme.hpp"
#include "../Constants.hpp"
#include "cplscheme/CouplingScheme.hpp"
#include "logging/LogMacros.hpp"

namespace precice::cplscheme::tests {

DummyCouplingScheme::DummyCouplingScheme(
    int numberIterations,
    int maxTimeSteps)
    : _numberIterations(numberIterations),
      _maxTimeSteps(maxTimeSteps)
{
}

void DummyCouplingScheme::initialize(
    double startTime,
    int    startTimesteps)
{
  PRECICE_ASSERT(not _isInitialized);
  _isInitialized = true;
  _isOngoing     = true;
  _timeWindows   = startTimesteps;
  _iterations    = 1;
}

CouplingScheme::ChangedMeshes DummyCouplingScheme::firstSynchronization(const CouplingScheme::ChangedMeshes &changes)
{
  PRECICE_ASSERT(_isInitialized);
  PRECICE_ASSERT(_isOngoing);
  PRECICE_ASSERT(changes.empty());
  return changes;
}

void DummyCouplingScheme::firstExchange()
{
  PRECICE_ASSERT(_isInitialized);
  PRECICE_ASSERT(_isOngoing);
}

CouplingScheme::ChangedMeshes DummyCouplingScheme::secondSynchronization()
{
  PRECICE_ASSERT(_isInitialized);
  PRECICE_ASSERT(_isOngoing);
  return {};
}

void DummyCouplingScheme::secondExchange()
{
  PRECICE_ASSERT(_isInitialized);
  PRECICE_ASSERT(_isOngoing);
  // Imagine we compute the convergence measure here
  _hasConverged = _iterations == _numberIterations;

  if (_hasConverged) {
    if (_timeWindows == _maxTimeSteps) {
      _isOngoing = false;
    }
    _timeWindows++;
    _iterations = 1;
  } else {
    _iterations++;
  }
  if (isImplicitCouplingScheme()) {
    PRECICE_DEBUG("advanced to {}-{}/{} (ongoing {})", _timeWindows, _iterations, _numberIterations, _isOngoing);
  } else {
    PRECICE_DEBUG("advanced to {} (ongoing {})", _timeWindows, _isOngoing);
  }
}

void DummyCouplingScheme::finalize()
{
  PRECICE_ASSERT(_isInitialized);
  PRECICE_ASSERT(not _isOngoing);
}

bool DummyCouplingScheme::isCouplingOngoing() const
{
  PRECICE_ASSERT(_isInitialized);
  if (_timeWindows <= _maxTimeSteps)
    return true;
  return false;
}

bool DummyCouplingScheme::isActionRequired(
    Action action) const
{
  if (!isImplicitCouplingScheme()) {
    PRECICE_DEBUG("return false (explicit)");
    return false;
  }
  if (action == CouplingScheme::Action::WriteCheckpoint) {
    if (_iterations == 1) {
      PRECICE_DEBUG("return true");
      return true;
    }
  }
  if (action == CouplingScheme::Action::ReadCheckpoint) {
    if (_iterations != 1) {
      PRECICE_DEBUG("return true");
      return true;
    }
  }
  PRECICE_DEBUG("return false");
  return false;
}

bool DummyCouplingScheme::hasConverged() const
{
  return _hasConverged;
}

} // namespace precice::cplscheme::tests
