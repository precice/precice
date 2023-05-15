#include "CouplingScheme.hpp"

namespace precice::cplscheme {

const double CouplingScheme::UNDEFINED_MAX_TIME = -1.0;

const int CouplingScheme::UNDEFINED_TIME_WINDOWS = -1;

const double CouplingScheme::UNDEFINED_TIME_WINDOW_SIZE = -1.0;

const int CouplingScheme::UNDEFINED_EXTRAPOLATION_ORDER = 0; // by default use zeroth order extrapolation.

const int CouplingScheme::UNDEFINED_MAX_ITERATIONS = -1;

std::string CouplingScheme::toString(Action action)
{
  static const std::map<CouplingScheme::Action, const char *> actionNames{
      {CouplingScheme::Action::WriteCheckpoint, "write-iteration-checkpoint"},
      {CouplingScheme::Action::ReadCheckpoint, "read-iteration-checkpoint"},
      {CouplingScheme::Action::InitializeData, "write-initial-data"}};

  return actionNames.at(action);
}

} // namespace precice::cplscheme
