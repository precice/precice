#include "utils/ManageUniqueIDs.hpp"

namespace precice
{
namespace utils
{

int ManageUniqueIDs::getFreeID()
{
  bool notFound = true;
  while (notFound) {
    if (_ids.count(_lowerLimit) == 0) {
      notFound = false;
    }
    _lowerLimit++;
  }
  _ids.insert(_lowerLimit - 1);
  return _lowerLimit - 1;
}

bool ManageUniqueIDs::insertID(int id)
{
  if (_ids.count(id) != 0)
    return false;
  else
    _ids.insert(id);
  return true;
}

void ManageUniqueIDs::resetIDs()
{
  _ids.clear();
  _lowerLimit = 0;
}

} // namespace utils
} // namespace precice
