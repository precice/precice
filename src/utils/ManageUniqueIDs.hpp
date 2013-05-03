// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_UTILS_MANAGEUNIQUEIDS_HPP_
#define PRECICE_UTILS_MANAGEUNIQUEIDS_HPP_

#include <set>

namespace precice {
namespace utils {

/**
 * @brief Manages a set of unique IDs.
 */
class ManageUniqueIDs
{
public:

   /**
    * @brief Constructor.
    */
   ManageUniqueIDs ();

   /**
    * @brief Returns the next free, i.e. unique, ID.
    */
   int getFreeID ();

   /**
    * @brief Inserts an ID which has to be unique.
    *
    * The inserted ID has to be different to all other IDs inserted and obtained
    * from getFreeID().
    */
   bool insertID ( int id );

   /**
    * @brief Resets all retrieved and inserted IDs.
    */
   void resetIDs ();

private:

   // @brief Stores all used IDs.
   std::set<int> _ids;

   // @brief Marks next ID to be given, from lower to higher values.
   int _lowerLimit;
};

}} // namespace precice, utils

#endif /* PRECICE_UTILS_MANAGEUNIQUEIDS_HPP_ */
