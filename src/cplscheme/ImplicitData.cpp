#include "cplscheme/ImplicitData.hpp"
#include "utils/assertion.hpp"

namespace precice::cplscheme {

void ImplicitData::add(DataID did, bool toKeep)
{
  PRECICE_ASSERT(entries.count(did) == 0);
  entries[did] = toKeep;
}

bool ImplicitData::contains(DataID did) const
{
  return entries.count(did) > 0;
}

bool ImplicitData::toKeep(DataID did) const
{
  PRECICE_ASSERT(contains(did));
  return entries.at(did);
}

} // namespace precice::cplscheme
