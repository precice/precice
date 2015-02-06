// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "PropertyContainer.hpp"
#include "utils/ManageUniqueIDs.hpp"

namespace precice {
namespace mesh {

const int PropertyContainer:: INDEX_GEOMETRY_ID = getFreePropertyID();

tarch::logging::Log PropertyContainer:: _log ( "precice::mesh::PropertyContainer" );

utils::ManageUniqueIDs * PropertyContainer:: _manageUniqueIDs = NULL;

PropertyContainer:: PropertyContainer ()
:
   _properties (),
   _parents ()
{}

const PropertyContainer & PropertyContainer:: getParent ( int index ) const
{
   assertion ( index >= 0 );
   assertion ( index < (int)_parents.size() );
   return *_parents[index];
}

bool PropertyContainer:: deleteProperty ( int properyID )
{
   if ( _properties.count(properyID) > 0 ) {
      _properties.erase (properyID);
      return true;
   }
   return false;
}


bool PropertyContainer:: hasProperty ( int propertyID ) const
{
   if ( _properties.find(propertyID) == _properties.end() ) {
      for ( size_t i=0; i < _parents.size(); i++ ) {
         if ( _parents[i]->hasProperty(propertyID) ) {
            return true;
         }

      }
      return false;
   }
   return true;
}

int PropertyContainer:: getFreePropertyID ()
{
   if ( _manageUniqueIDs == NULL ) {
      _manageUniqueIDs = new utils::ManageUniqueIDs ();
   }
   return _manageUniqueIDs->getFreeID();
}

void PropertyContainer:: resetPropertyIDCounter ()
{
   if ( _manageUniqueIDs != NULL ) {
      _manageUniqueIDs->resetIDs ();
   }
}

}}
