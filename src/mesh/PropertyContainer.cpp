#include "PropertyContainer.hpp"
#include "utils/ManageUniqueIDs.hpp"

namespace precice
{
namespace mesh
{

const int PropertyContainer::INDEX_GEOMETRY_ID = getFreePropertyID();

logging::Logger PropertyContainer::_log("precice::mesh::PropertyContainer");

utils::ManageUniqueIDs *PropertyContainer::_manageUniqueIDs = nullptr;

PropertyContainer::PropertyContainer()
    : _properties(),
      _parents()
{
}

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
  if (_manageUniqueIDs == nullptr) {
    _manageUniqueIDs = new utils::ManageUniqueIDs();
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
