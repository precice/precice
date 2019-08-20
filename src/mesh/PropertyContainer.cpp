#include "PropertyContainer.hpp"
#include "utils/ManageUniqueIDs.hpp"

namespace precice
{
namespace mesh
{

const int PropertyContainer::INDEX_GEOMETRY_ID = getFreePropertyID();

std::unique_ptr<utils::ManageUniqueIDs> PropertyContainer::_manageUniqueIDs;
logging::Logger PropertyContainer::_log{"mesh::PropertyContainer"};

const PropertyContainer &PropertyContainer::getParent(size_t index) const
{
  return *_parents.at(index);
}

bool PropertyContainer::deleteProperty(int propertyID)
{
  if (_properties.count(propertyID) > 0) {
    _properties.erase(propertyID);
    return true;
  }
  return false;
}

bool PropertyContainer::hasProperty(int propertyID) const
{
  if (_properties.find(propertyID) == _properties.end()) {
    for (const auto &elem : _parents) {
      if (elem->hasProperty(propertyID)) {
        return true;
      }
    }
    return false;
  }
  return true;
}

int PropertyContainer::getFreePropertyID()
{
  if (not _manageUniqueIDs) {
    _manageUniqueIDs.reset(new utils::ManageUniqueIDs);
  }
  return _manageUniqueIDs->getFreeID();
}

void PropertyContainer::resetPropertyIDCounter()
{
  if (_manageUniqueIDs != nullptr) {
    _manageUniqueIDs->resetIDs();
  }
}
}
}
