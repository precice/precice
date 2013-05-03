// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "utils/ManageUniqueIDs.hpp"

namespace precice {
namespace utils {

ManageUniqueIDs:: ManageUniqueIDs ()
:
   _ids (),
   _lowerLimit (0)
{}

int ManageUniqueIDs:: getFreeID ()
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

bool ManageUniqueIDs:: insertID ( int id )
{
   if (_ids.count(id) != 0)
      return false;
   else
      _ids.insert (id);
   return true;
}

void ManageUniqueIDs:: resetIDs ()
{
   _ids.clear ();
   _lowerLimit = 0;
}

}} // namespace precice, utils
