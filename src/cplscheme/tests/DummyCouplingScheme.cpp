#include "DummyCouplingScheme.hpp"
#include "../Constants.hpp"
#include "cplscheme/CouplingScheme.hpp"
#include "logging/LogMacros.hpp"

namespace precice::cplscheme::tests {

DummyCouplingScheme::DummyCouplingScheme(
    int numberIterations,
    int maxTimesteps)
    : _numberIterations(numberIterations),
      _maxTimesteps(maxTimesteps)
{
}

void DummyCouplingScheme::initialize(
    double startTime,
    int    startTimesteps)
{
  PRECICE_ASSERT(not _isInitialized);
  _isInitialized = true;
  _isOngoing     = true;
  _timesteps     = startTimesteps;
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

bool DummyCouplingScheme::secondExchange()
{
  PRECICE_ASSERT(_isInitialized);
  PRECICE_ASSERT(_isOngoing);
  if (_iterations == _numberIterations) {
    if (_timesteps == _maxTimesteps) {
      _isOngoing = false;
    }
    _timesteps++;
    _iterations = 0;
  }
  _iterations++;
  return true;
}

void DummyCouplingScheme::finalize()
{
  PRECICE_ASSERT(_isInitialized);
  PRECICE_ASSERT(not _isOngoing);
}

bool DummyCouplingScheme::isCouplingOngoing() const
{
  PRECICE_ASSERT(_isInitialized);
  if (_timesteps <= _maxTimesteps)
    return true;
  return false;
}

bool DummyCouplingScheme::isActionRequired(
    const std::string &actionName) const
{
  if (_numberIterations > 1) {
    if (actionName == constants::actionWriteIterationCheckpoint()) {
      if (_iterations == 1) {
        PRECICE_DEBUG("return true");
        return true;
      }
    } else if (actionName == constants::actionReadIterationCheckpoint()) {
      if (_iterations != 1) {
        PRECICE_DEBUG("return true");
        return true;
      }
    }
  }
  PRECICE_DEBUG("return false");
  return false;
}

} // namespace precice::cplscheme::tests
