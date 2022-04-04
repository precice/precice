#include <acceleration/Acceleration.hpp>
#include "cplscheme/CouplingData.hpp"
#include "utils/Helpers.hpp"

namespace precice {
namespace acceleration {

void Acceleration::checkDataIDs(const DataMap &cplData) const
{
#ifndef NDEBUG
  for (int id : getDataIDs()) {
    bool valid = utils::contained(id, cplData);
    PRECICE_ASSERT(valid, "Data with ID {} unknown.", id);
  }
#endif
}

} // namespace acceleration
} // namespace precice
