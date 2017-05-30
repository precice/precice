
#include "DummyCouplingScheme.hpp"
#include "../Constants.hpp"


namespace precice {
namespace cplscheme {
namespace tests {

logging::Logger DummyCouplingScheme::
   _log("cplscheme::tests::DummyCouplingScheme");

DummyCouplingScheme:: DummyCouplingScheme
(
  int numberIterations,
  int maxTimesteps )
:
  _numberIterations(numberIterations),
  _iterations(0),
  _maxTimesteps(maxTimesteps),
  _timesteps(0),
  _isInitialized(false),
  _isOngoing(false)
{}

void DummyCouplingScheme:: initialize
(
  double startTime,
  int    startTimesteps )
{
  preciceTrace("initialize()", startTime, startTimesteps);
  assertion(not _isInitialized);
  _isInitialized = true;
  _isOngoing = true;
  _timesteps = startTimesteps;
  _iterations=1;
}

void DummyCouplingScheme:: advance()
{
  preciceTrace("advance()", _iterations, _timesteps);
  assertion(_isInitialized);
  assertion(_isOngoing);
  if (_iterations == _numberIterations){
    if (_timesteps == _maxTimesteps){
      _isOngoing = false;
    }
    _timesteps++;
    _iterations = 0;
  }
  _iterations++;
}

void DummyCouplingScheme:: finalize()
{
  preciceTrace("finalize()");
  assertion(_isInitialized);
  assertion(not _isOngoing);
}

bool DummyCouplingScheme:: isCouplingOngoing() const
{
  if (_timesteps <= _maxTimesteps) return true;
  return false;
}

bool DummyCouplingScheme:: isActionRequired
(
  const std::string& actionName ) const
{
  preciceTrace("isActionRequired()", actionName);
  if (_numberIterations > 1){
    if (actionName == constants::actionWriteIterationCheckpoint()){
      if (_iterations == 1) {
        DEBUG("return true");
        return true;
      }
    }
    else if (actionName == constants::actionReadIterationCheckpoint()){
      if (_iterations != 1) {
        DEBUG("return true");
        return true;
      }
    }
  }
  DEBUG("return false");
  return false;
}

}}} // namespace precice, cplscheme, tests
