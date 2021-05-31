#include "utils/ManageUniqueIDs.hpp"
#include <boost/container/detail/flat_tree.hpp>

namespace precice {
namespace utils {

int ManageUniqueIDs::getFreeID()
{
  bool notFound = true;
  while (notFound) {
    if (_ids.find(_lowerLimit) == _ids.end()) {
      notFound = false;
    }
    _lowerLimit++;
  }
  _ids.insert(_ids.end(), _lowerLimit - 1);
  return _lowerLimit - 1;
}

bool ManageUniqueIDs::insertID(int id)
{
  if (_ids.find(id) != _ids.end()) {
    return false;
  }

  _ids.insert(_ids.end(), id);
  return true;
}

void ManageUniqueIDs::resetIDs()
{
  _ids.clear();
  _lowerLimit = 0;
}

} // namespace utils
} // namespace precice
