// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_LA_TRAITS_EQUALSSCALARS_H_
#define _TARCH_LA_TRAITS_EQUALSSCALARS_H_

#include "tarch/la/traits/IsVector.h"
#include "tarch/la/traits/VectorTraits.h"
#include "tarch/la/traits/IsMatrix.h"
#include "tarch/la/traits/MatrixTraits.h"

namespace tarch {
namespace la {

//template<typename VectorA, typename VectorB,
//         bool Decider = IsVector<VectorA>::value && IsVector<VectorB>::value>
template<typename VectorA, typename VectorB, typename Enabler = void>
struct EqualScalars
{
  static const int value = false;
};

template<typename VectorA, typename VectorB>
struct EqualScalars<
  VectorA, VectorB,
  typename std::enable_if<IsVector<VectorA>::value && IsVector<VectorB>::value>::type>
{
  typedef typename VectorTraits<VectorA>::Scalar ScalarA;
  typedef typename VectorTraits<VectorB>::Scalar ScalarB;
  static const int value = std::is_same<ScalarA,ScalarB>::value;
};

template<typename MatrixA, typename MatrixB>
struct EqualScalars<
  MatrixA, MatrixB,
  typename std::enable_if<IsMatrix<MatrixA>::value && IsMatrix<MatrixB>::value>::type>
{
  typedef typename MatrixTraits<MatrixA>::Scalar ScalarA;
  typedef typename MatrixTraits<MatrixB>::Scalar ScalarB;
  static const int value = std::is_same<ScalarA,ScalarB>::value;
};

}} // namespace tarch, la

#endif /* _TARCH_LA_TRAITS_EQUALSSCALARS_H_ */
