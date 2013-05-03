// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_UTILS_POINTERVECTOR_HPP_
#define PRECICE_UTILS_POINTERVECTOR_HPP_

#include "Globals.hpp"
#include <vector>
#include "boost/iterator/indirect_iterator.hpp"

namespace precice {
namespace utils {

/**
 * @brief Wrapper around std::vector for transparent handling of pointer objects.
 */
template< typename CONTENT_T >
class ptr_vector
{
public:

   // @brief Type of wrapped vector container.
   typedef std::vector<CONTENT_T*> container;
   typedef CONTENT_T value_type; // necessary to be standard conform

   // @brief Type of iterator hiding pointers.
   typedef boost::indirect_iterator<
              typename container::iterator,
              CONTENT_T*,
              boost::use_default,
              CONTENT_T&
           > iterator;

   // @brief Type of const_iterator hiding pointers.
   typedef boost::indirect_iterator<
              typename container::const_iterator,
              CONTENT_T *,
              boost::use_default,
              CONTENT_T &
           > const_iterator;

   /**
    * @brief Returns iterator to first element in vector.
    */
   iterator begin ()
   {
      return iterator(_content.begin());
   }

   /**
    * @brief Returns const interator to first element in vector.
    */
   const_iterator begin () const
   {
      return const_iterator(_content.begin());
   }

   /**
    * @brief Returns iterator behind last element in vector.
    */
   iterator end ()
   {
      return iterator(_content.end());
   }

   /**
    * @brief Returns const_iterator behind last element in vector.
    */
   const_iterator end () const
   {
      return const_iterator(_content.end());
   }

   /**
    * @brief Returns count of elements in vector.
    */
   size_t size () const
   {
      return _content.size ();
   }

   /**
    * @brief Returns reference to element with given index [0, count[.
    */
   CONTENT_T & operator[] ( size_t index )
   {
      assertion ( index >= 0 );
      assertion ( index < _content.size() );
      return *_content[index];
   }

   /**
    * @brief Returns const reference to element with given index [0, count[.
    */
   const CONTENT_T & operator[] ( size_t index ) const
   {
      assertion ( index >= 0 );
      assertion ( index < _content.size() );
      return *_content[index];
   }

   CONTENT_T& back()
   {
     return *_content.back();
   }

   const CONTENT_T& back() const
   {
     return *_content.back();
   }

   /**
    * @brief Adds element to the end of the vector.
    */
   void push_back ( CONTENT_T * content )
   {
      _content.push_back ( content );
   }

   /**
    * @brief Inserts elements into vector.
    *
    * @param position [IN] Iterator at position before which the elements are
    *        inserted.
    * @param from [IN] Iterator at position of first element to be inserted.
    * @param to [IN] Iterator at position behind last element to be inserted.
    */
   void insert ( iterator position, iterator from, iterator to  )
   {
      _content.insert ( position.base(), from.base(), to.base() );
   }

   /**
    * @brief Removes all pointers to elements from the vector (no deletetion).
    */
   void clear ()
   {
      _content.clear ();
   }

   /**
    * @brief Returns true, if no pointers are contained in the vector.
    */
   bool empty () const
   {
      return _content.empty ();
   }

   void deleteElements ()
   {
      foreach ( CONTENT_T * elem, _content ) {
         assertion ( elem != NULL );
         delete ( elem );
      }
   }

private:

   container _content;
};

// --------------------------------------------------------- HEADER DEFINITIONS

}} // namespace precice, utils

#endif /* PRECICE_UTILS_POINTERVECTOR_HPP_ */
