// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_LA_TRAITS_ISVECTOR_H_
#define _TARCH_LA_TRAITS_ISVECTOR_H_

#include <vector>
#include <deque>

namespace tarch {
   namespace la {
      template<int Size, typename Scalar> class Vector;
      template<typename Scalar> class DynamicVector;
      template<int Size, typename Scalar> class WrappedVector;
   }
}


namespace tarch {
namespace la {

template< typename VECTOR_T >
struct IsVector
{
   static const int value = false;
};

template<int Size, typename Scalar>
struct IsVector<Vector<Size,Scalar> >
{
   static const int value = true;
};

template<typename Scalar>
struct IsVector<DynamicVector<Scalar> >
{
   static const int value = true;
};

template<int Size, typename Scalar>
struct IsVector<WrappedVector<Size,Scalar> >
{
   static const int value = true;
};

template< typename Scalar >
struct IsVector<std::vector<Scalar> >
{
   static const int value = true;
};

template< typename Scalar >
struct IsVector<std::deque<Scalar> >
{
   static const int value = true;
};

}} // namespace tarch, la

#endif /* _TARCH_LA_TRAITS_ISVECTOR_H_ */
