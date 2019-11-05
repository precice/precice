#pragma once

#include "utils/assertion.hpp"
#include <boost/any.hpp>
#include <map>
#include <vector>

namespace precice
{
namespace utils
{
class ManageUniqueIDs;
}
}

namespace precice
{
namespace mesh
{

/**
 * @brief Implements hierarchical dynamical properties of arbitrary type.
 *
 * This class inherited by other classes to make them have properties. By using
 * the boost::any type, the type of a property can be arbitrary, but access is
 * typed at compile time. A property is accessed via it's ID, which has to be
 * unique among all property IDs of one PropertyContainer object. Properties
 * are created and deleted dynamically. Hierarchical behavior is introduced by
 * parent pointers to higher level PropertyContainers. There can be multiple
 * parents.
 */
class PropertyContainer
{
public:
  PropertyContainer& operator=(PropertyContainer &&) = delete;

  virtual ~PropertyContainer(){};

  // Shortform for the type of a property.
  using PropertyType = boost::any;

  /// ID for the property labeling geometry IDs.
  static const int INDEX_GEOMETRY_ID;

  /// Returns a new, globally unique property ID.
  static int getFreePropertyID();

  /**
    * @brief Resets the ID counter of all properties.
    *
    * This function is meant for test cases. Already created PropertyContainer
    * objects might be in inconsistent state after calling this function.
    */
  static void resetPropertyIDCounter();

  /// Enables hierarchical property behavior.
  void addParent(PropertyContainer &parent)
  {
    _parents.push_back(&parent);
  }

  /// Returns the number of parents.
  int getParentCount() const
  {
    return static_cast<int>(_parents.size());
  }

  /// Returns the parent corresponding to the given index (0 ... count).
  const PropertyContainer &getParent(size_t index) const;

  /**
     * @brief Create (if not existing) and set property value.
     *
     * @param[in] propertyID ID of the property.
     * @param[in] value Value to be set for the property.
     */
  template <typename value_t>
  void setProperty(int propertyID, const value_t &value)
  {
    _properties[propertyID] = value;
  }

  /**
     * @brief Returns true, when a property with given ID can be retrieved.
     *
     * True is also returned, when one or more parent PropertyContainer are set
     * and have the property with given ID (once or more).
     */
  bool hasProperty(int propertyID) const;

  /**
     * @brief Deletes a set property.
     *
     * This allows temporary variables for container objects.
     * Gives an error, when the property does not exist.
     *
     * @param[in] propertyID ID of the property to delete
     * @return true, if property existed, false if not
     */
  bool deleteProperty(int propertyID);

  /**
   * @brief Returns the value of the property with given ID.
   *
   * @pre The property has to exist for the object, or for its PropertyContainer object
   * @pre The type of the property has to coincide with the one specified as
   *   explicit template parameter when calling getProperty<type>().
   */
  template <typename value_t>
  const value_t &getProperty(int propertyID) const;
  
  /**
   * @brief Recursively looks up a property ID.
   *
   * @param[in] propertyID ID to lookup
   * @param[out] properties found
   */
  template <typename value_t>
  void getProperties(int propertyID, std::vector<value_t> &properties) const;

private:
  /** The logger for the PropertyContainer.
   *
   * @attention There will be many instances of this class.
   * Thus using a static logger is crucial, as creating a logger is a costly
   * and an allocation-heavy operation,
   */
  static logging::Logger _log;

  /// Manager to ensure unique identification of all properties.
  static std::unique_ptr<utils::ManageUniqueIDs> _manageUniqueIDs;

  /// Properties (local for every instance)
  std::map<int, PropertyType> _properties;

  /// Must be set, if hierarchical properties are wanted
  std::vector<PropertyContainer *> _parents;
};

// --------------------------------------------------------- HEADER DEFINITIONS

template <typename value_t>
const value_t &PropertyContainer::getProperty(int propertyID) const
{
  std::map<int, PropertyType>::const_iterator iter;
  iter = _properties.find(propertyID);
  if (iter == _properties.end()) {
    for (const auto &prop : _parents) {
      if (prop->hasProperty(propertyID)) {
        return prop->getProperty<value_t>(propertyID);
      }
    }
    PRECICE_ERROR("No property with id = " << propertyID);
  }
  PRECICE_ASSERT(not iter->second.empty());
  // When the type of value_t does not match that of the any, NULL is returned.
  PRECICE_ASSERT(boost::any_cast<value_t>(&iter->second) != nullptr);
  return *boost::any_cast<value_t>(&iter->second);
}

template <typename value_t>
void PropertyContainer::getProperties(int propertyID, std::vector<value_t> &properties) const
{
  auto iter = _properties.find(propertyID);
  if (iter != _properties.end()) {
    PRECICE_ASSERT(not iter->second.empty());
    // When the type of value_t does not match that of the any, NULL is returned.
    PRECICE_ASSERT(boost::any_cast<value_t>(&iter->second) != nullptr);
    properties.push_back(boost::any_cast<value_t>(iter->second));
  } else {
    for (size_t i = 0; i < _parents.size(); i++) {
      _parents[i]->getProperties(propertyID, properties);
    }
  }
}
}
} // namespace precice, mesh
