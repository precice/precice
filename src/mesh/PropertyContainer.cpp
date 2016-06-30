#include "PropertyContainer.hpp"
#include "utils/ManageUniqueIDs.hpp"

namespace precice {
namespace mesh {

const int PropertyContainer:: INDEX_GEOMETRY_ID = getFreePropertyID();

tarch::logging::Log PropertyContainer:: _log ( "precice::mesh::PropertyContainer" );

utils::ManageUniqueIDs * PropertyContainer:: _manageUniqueIDs = nullptr;

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
      for (auto & elem : _parents) {
         if ( elem->hasProperty(propertyID) ) {
            return true;
         }

      }
      return false;
   }
   return true;
}

int PropertyContainer:: getFreePropertyID ()
{
   if ( _manageUniqueIDs == nullptr ) {
      _manageUniqueIDs = new utils::ManageUniqueIDs ();
   }
   return _manageUniqueIDs->getFreeID();
}

void PropertyContainer:: resetPropertyIDCounter ()
{
   if ( _manageUniqueIDs != nullptr ) {
      _manageUniqueIDs->resetIDs ();
   }
}

}}
