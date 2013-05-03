// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_CONTAINER_PROPERTYCONTAINER_HPP_
#define PRECICE_CONTAINER_PROPERTYCONTAINER_HPP_

#include "utils/Globals.hpp"
#include "utils/Dimensions.hpp"
#include "boost/any.hpp"
#include <map>
#include <vector>

namespace precice {
   namespace utils {
      class ManageUniqueIDs;
   }
}

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace mesh {



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

   /**
    * @brief Shortform for the type of a property.
    */
   typedef boost::any PropertyType;

   // @brief ID for the property labeling geometry IDs.
   static const int INDEX_GEOMETRY_ID;

   /**
    * @brief Returns a new, globally unique property ID.
    */
   static int getFreePropertyID();

   /**
    * @brief Resets the ID counter of all properties.
    *
    * This function is meant for test cases. Already created PropertyContainer
    * objects might be in inconsistent state after calling this function.
    */
   static void resetPropertyIDCounter();

   /**
    * @brief Constructor.
    */
   PropertyContainer();

   /**
    * @brief Destructor.
    */
   virtual ~PropertyContainer() {};

    /**
     * @brief Enables hierarchical property behavior.
     */
    void addParent ( PropertyContainer& parent )
    {
       _parents.push_back(& parent);
    }

    /**
     * @brief Returns the number of parents.
     */
    int getParentCount() const
    {
       return (int)_parents.size();
    }

    /**
     * @brief Returns the parent corresponding to the given index (0 ... count).
     */
    const PropertyContainer& getParent ( int index ) const;

    /**
     * @brief Create (if not existing) and set property value.
     *
     * @param propertyID [IN] ID of the property.
     * @param value [IN] Value to be set for the property.
     */
    template< typename value_t >
    void setProperty ( int             propertyID,
                       const value_t & value )
    {
       _properties[propertyID] = value;
    }

    /**
     * @brief Returns true, when a property with given ID can be retrieved.
     *
     * True is also returned, when one or more parent PropertyContainer are set
     * and have the property with given ID (once or more).
     */
    bool hasProperty ( int propertyID ) const;

    /**
     * @brief Deletes a set property.
     *
     * This allows temporary variables for container objects.
     * Gives an error, when the property does not exist.
     *
     * @param index [IN] Index of the property to delete
     * @return true, if property existed, false if not
     */
    bool deleteProperty ( int propertyID );

    /**
     * @brief Returns the value of the property with given ID.
     *
     * Prerequesits:
     * - the property has to exist for the object, or for its parent
     *   PropertyContainer object
     * - the type of the property has to coincide with the one specified as
     *   explicit template parameter when calling getProperty<type>().
     */
    template< typename value_t >
    const value_t& getProperty ( int propertyID ) const;

    /**
     * @brief Returns all properties of this and parent PropertyContainer objects.
     */
    template< typename value_t>
    void getProperties ( int                    propertyID,
                         std::vector<value_t> & properties );

private:

   // @brief Logging device.
   static tarch::logging::Log _log;

   // @brief Manager to ensure unique identification of all properties.
   static utils::ManageUniqueIDs * _manageUniqueIDs;

   // @brief Properties (local for every instance)
   std::map<int, PropertyType> _properties;

   // @brief Must be set, if hierarchical properties are wanted
   std::vector<PropertyContainer *> _parents;
};



// --------------------------------------------------------- HEADER DEFINITIONS



template< typename value_t >
const value_t& PropertyContainer:: getProperty
(
  int propertyID ) const
{
  std::map<int, PropertyType>::const_iterator iter;
  iter = _properties.find ( propertyID );
  if ( iter == _properties.end() ) {
    for ( size_t i=0; i < _parents.size(); i++ ) {
      if ( _parents[i]->hasProperty(propertyID) ) {
        return _parents[i]->getProperty<value_t> ( propertyID );
      }
    }
    preciceError ( "getProperty()", "No property with id = " << propertyID );
  }
  assertion ( not iter->second.empty() );
  // When the type of value_t does not match that of the any, NULL is returned.
  assertion ( boost::any_cast<value_t>(&iter->second) != NULL );
  return * boost::any_cast<value_t> ( & iter->second );
}

template< typename value_t >
void PropertyContainer:: getProperties
(
  int                   propertyID,
  std::vector<value_t>& properties )
{
  std::map<int, PropertyType>::const_iterator iter;
  iter = _properties.find(propertyID);
  if (iter != _properties.end()){
    assertion(not iter->second.empty());
    // When the type of value_t does not match that of the any, NULL is returned.
    assertion(boost::any_cast<value_t>(&iter->second) != NULL);
    properties.push_back(boost::any_cast<value_t>(iter->second));
  }
  else {
    for (size_t i=0; i < _parents.size(); i++){
      _parents[i]->getProperties(propertyID, properties);
    }
  }
}

}} // namespace precice, mesh

#endif /* PRECICE_CONTAINER_PROPERTYCONTAINER_HPP_ */
