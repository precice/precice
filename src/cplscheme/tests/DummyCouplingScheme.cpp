// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License

#include "DummyCouplingScheme.hpp"
#include "../Constants.hpp"


namespace precice {
namespace cplscheme {
namespace tests {

tarch::logging::Log DummyCouplingScheme::
   _log("precice::cplscheme::tests::DummyCouplingScheme");

DummyCouplingScheme:: DummyCouplingScheme
(
  int numberIterations,
  int maxTimesteps )
:
  _numberIterations(numberIterations),
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
  preciceTrace2("initialize()", startTime, startTimesteps);
  assertion(not _isInitialized);
  _isInitialized = true;
  _isOngoing = true;
  _timesteps = startTimesteps;
}

void DummyCouplingScheme:: advance()
{
  preciceTrace("advance()");
  assertion(_isInitialized);
  assertion(_isOngoing);
  _iterations++;
  if (_iterations == _numberIterations){
    _timesteps++;
    if (_timesteps == _maxTimesteps){
      _isOngoing = false;
    }
    _iterations = 0;
  }
}

void DummyCouplingScheme:: finalize()
{
  preciceTrace("finalize()");
  assertion(_isInitialized);
  assertion(not _isOngoing);
}

bool DummyCouplingScheme:: isActionRequired
(
  const std::string& actionName ) const
{
  preciceTrace1("isActionRequired()", actionName);
  if (_numberIterations > 1){
    if (actionName == constants::actionWriteIterationCheckpoint()){
      if (_iterations == 0) {
        return true;
      }
    }
    else if (actionName == constants::actionReadIterationCheckpoint()){
      if (_iterations < _numberIterations) {
        return true;
      }
    }
  }
  return false;
}

}}} // namespace precice, cplscheme, tests
