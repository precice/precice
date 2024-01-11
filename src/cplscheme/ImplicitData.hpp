#pragma once

#include <map>

#include "precice/impl/Types.hpp"

namespace precice::cplscheme {

struct ImplicitData {
  void add(DataID did, bool toKeep);

  bool contains(DataID did) const;
  bool toKeep(DataID did) const;

  std::map<DataID, bool> entries;
};

} // namespace precice::cplscheme
