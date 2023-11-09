#include "CouplingScheme.hpp"

#include <algorithm>

namespace precice::cplscheme {

const double CouplingScheme::UNDEFINED_MAX_TIME = -1.0;

const int CouplingScheme::UNDEFINED_TIME_WINDOWS = -1;

const double CouplingScheme::UNDEFINED_TIME_WINDOW_SIZE = -1.0;

const int CouplingScheme::UNDEFINED_MAX_ITERATIONS = -1;
const int CouplingScheme::INFINITE_MAX_ITERATIONS  = -2;
const int CouplingScheme::UNDEFINED_MIN_ITERATIONS = -1;

std::string CouplingScheme::toString(Action action)
{
  static const std::map<CouplingScheme::Action, const char *> actionNames{
      {CouplingScheme::Action::WriteCheckpoint, "write-iteration-checkpoint"},
      {CouplingScheme::Action::ReadCheckpoint, "read-iteration-checkpoint"},
      {CouplingScheme::Action::InitializeData, "write-initial-data"}};

  return actionNames.at(action);
}

CouplingScheme::ExchangePlan &CouplingScheme::ExchangePlan::tidy()
{
  auto sortUnique = [](auto &c) {
    std::sort(c.begin(), c.end());
    c.erase(std::unique(c.begin(), c.end()), c.end());
  };
  sortUnique(sendExplicit);
  sortUnique(sendImplicit);
  sortUnique(receiveExplicit);
  sortUnique(receiveImplicit);
  return *this;
}

CouplingScheme::ExchangePlan &CouplingScheme::ExchangePlan::operator+=(const ExchangePlan &other)
{
  auto append = [](auto &dst, const auto &src) { dst.insert(dst.end(), src.begin(), src.end()); };
  append(sendExplicit, other.sendExplicit);
  append(sendImplicit, other.sendImplicit);
  append(receiveExplicit, other.receiveExplicit);
  append(receiveImplicit, other.receiveImplicit);
  return tidy();
}

} // namespace precice::cplscheme
