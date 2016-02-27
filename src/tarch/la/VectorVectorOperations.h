// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_LA_VECTORVECTOROPERATIONS_H_
#define _TARCH_LA_VECTORVECTOROPERATIONS_H_

#include "tarch/la/traits/IsVector.h"
#include "tarch/la/traits/VectorTraits.h"
#include "tarch/la/traits/DeduceScalar.h"
//#include "tarch/la/traits/EqualScalars.h"
#include "tarch/la/VectorAssign.h"
#include "tarch/utils/EnableIf.h"

namespace tarch {
  namespace la {

    /**
     * Adds rVector to lVector and assigns the result to lVector.
     *
     * No temporary vector is created.
     */
    template<typename LVector, typename RVector>
      typename std::enable_if<IsVector<LVector>::value && IsVector<RVector>::value /*&& EqualScalars<LVector,RVector>::value*/,
      LVector&
    >::type operator+= (
      LVector&       lVector,
      const RVector& rVector
    );

    /**
     * Subtracts rVector from lVector and assigns the result to lVector.
     *
     * No temporary vector is created.
     */
    template<typename LVector, typename RVector>
      typename std::enable_if<IsVector<LVector>::value && IsVector<RVector>::value /*&& EqualScalars<LVector,RVector>::value*/,
      LVector&
    >::type operator-= (
      LVector&       lVector,
      const RVector& rVector
    );

    /**
     * Adds rVector to lVector.
     *
     * A temporary vector is created and copied to store return back the result.
     */
    template<typename LVector, typename RVector>
      typename std::enable_if<
      IsVector<LVector>::value && IsVector<RVector>::value /*&& EqualScalars<LVector,RVector>::value*/,
      LVector
    >::type operator+ (
      const LVector& lVector,
      const RVector& rVector
    );

    /**
     * Subtracts rVector from lVector.
     *
     * A temporary vector is created and copied to store return back the result.
     */
    template<typename LVector, typename RVector>
      typename std::enable_if<
      IsVector<LVector>::value && IsVector<RVector>::value /*&& EqualScalars<LVector,RVector>::value*/,
      LVector
    >::type operator- (
      const LVector& lVector,
      const RVector& rVector
    );

    /**
     * Multiplies every component of the vectors with each other and writes the
     * results into result.
     */
    template<typename LVector, typename RVector, typename ResultVector>
      typename std::enable_if<
      IsVector<LVector>::value && IsVector<RVector>::value && IsVector<ResultVector>::value /*&& EqualScalars<LVector,RVector>::value*/,
      ResultVector&
    >::type multiplyComponents (
      const LVector& lVector,
      const RVector& rVector,
      ResultVector&  result
    );

    /**
     * Performs the dot (=inner) product of two vectors.
     */
    template<typename LVector, typename RVector>
      typename utils::LazyEnableIf<
      IsVector<LVector>::value && IsVector<RVector>::value /*&& EqualScalars<LVector,RVector>::value*/,
      utils::LazyType<typename VectorTraits<LVector>::Scalar>
    >::Type operator* (
      const LVector& lVector,
      const RVector& rVector
    );

    /**
     * Divides every component of the vectors with each other and returns
     * the resulting vector.
     */
    template<typename LVector, typename RVector>
    typename std::enable_if<
      IsVector<LVector>::value && IsVector<RVector>::value,
      LVector
    >::type operator/ (
      const LVector& lVector,
      const RVector& rVector
    );

    /**
     * Performs the dot (=inner) product of two vectors.
     */
    template<typename LVector, typename RVector>
      typename utils::LazyEnableIf<
      IsVector<LVector>::value && IsVector<RVector>::value /*&& EqualScalars<LVector,RVector>::value*/,
      utils::LazyType<typename VectorTraits<LVector>::Scalar>
    >::Type dot (
      const LVector & lVector,
      const RVector & rVector
    );

    /**
     * Performs the cross product of two 3D vectors into result.
     */
    template<typename Vector>
      typename std::enable_if< IsVector<Vector>::value,
      Vector&
    >::type cross (
      const Vector& lVector,
      const Vector& rVector,
      Vector&       result
    );

    /**
     * Compares to vectors on equality by means of a numerical accuracy.
     */
    template<typename LVector, typename RVector>
      typename std::enable_if<
      IsVector<RVector>::value && IsVector<LVector>::value /*&& EqualScalars<LVector,RVector>::value*/,
      bool
    >::type equals (
      const LVector&                         lVector,
      const RVector&                         rVector,
      typename VectorTraits<LVector>::Scalar tolerance = NUMERICAL_ZERO_DIFFERENCE
    );

    /**
     * Compares sequentially every component pair of lVector and rVector, and stops
     * as soon as one is greater than the other.
     *
     * Defines an absolute pairwise ordering between (unequal) vectors.
     */
    template<typename LVector, typename RVector>
      typename std::enable_if<
      IsVector<RVector>::value && IsVector<LVector>::value /*&& EqualScalars<RVector,LVector>::value*/,
      bool
    >::type firstGreater (
      const LVector&                         lVector,
      const RVector&                         rVector,
      typename VectorTraits<LVector>::Scalar tolerance = NUMERICAL_ZERO_DIFFERENCE
    );

    /**
     * Returns true, if one component of lVector is greater than the corresponding
     * component in rVector up to numerical accuracy.
     */
    template<typename LVector, typename RVector>
      typename std::enable_if<
      IsVector<RVector>::value && IsVector<LVector>::value /*&& EqualScalars<RVector,LVector>::value*/,
      bool
    >::type oneGreater (
      const LVector&                         lVector,
      const RVector&                         rVector,
      typename VectorTraits<LVector>::Scalar tolerance = NUMERICAL_ZERO_DIFFERENCE
    );

    /**
     * Returns true, if one component of lVector is greater or equals than the
     * corresponding component in rVector up to numerical accuracy.
     */
    template<typename LVector, typename RVector>
      typename std::enable_if<
      IsVector<RVector>::value && IsVector<LVector>::value /*&& EqualScalars<RVector,LVector>::value*/,
      bool
    >::type oneGreaterEquals (
      const LVector&                         lVector,
      const RVector&                         rVector,
      typename VectorTraits<LVector>::Scalar tolerance = NUMERICAL_ZERO_DIFFERENCE
    );

    /**
     * Returns true, if all components of lVector are greater than the corresponding
     * components in rVector up to numerical accuracy.
     */
    template<typename LVector, typename RVector>
      typename std::enable_if<
      IsVector<RVector>::value && IsVector<LVector>::value /*&& EqualScalars<RVector,LVector>::value*/,
      bool
    >::type allGreater (
      const LVector&                         lVector,
      const RVector&                         rVector,
      typename VectorTraits<LVector>::Scalar tolerance = NUMERICAL_ZERO_DIFFERENCE
    );

    /**
     * Returns true, if all components of lVector are greater or equals than the corresponding
     * components in rVector up to numerical accuracy.
     */
    template<typename LVector, typename RVector>
      typename std::enable_if<
      IsVector<RVector>::value && IsVector<LVector>::value,
      bool
    >::type allGreaterEquals (
      const LVector&                         lVector,
      const RVector&                         rVector,
      typename VectorTraits<LVector>::Scalar tolerance = NUMERICAL_ZERO_DIFFERENCE
    );

    /**
     * Bit-wise comparison of the components of two vectors for equality.
     *
     * This method should not be used for floating-point valued vectors. Instead,
     * equals() is the suitable comparison.
     */
    template<typename LVector, typename RVector>
      typename std::enable_if<
      IsVector<RVector>::value && IsVector<LVector>::value /*&& EqualScalars<LVector,RVector>::value*/,
      bool
    >::type operator== (
      const LVector& lVector,
      const RVector& rVector
    );

    /**
     * Bit-wise comparison of the components of two vectors for inequality.
     *
     * This method should not be used for floating-point valued vectors. Instead,
     * !equals() is the suitable comparison.
     */
    template<typename LVector, typename RVector>
      typename std::enable_if<
      IsVector<RVector>::value && IsVector<LVector>::value /*&& EqualScalars<LVector,RVector>::value*/,
      bool
    >::type operator!= (
      const LVector& lVector,
      const RVector& rVector
    );

    /**
     * Return Index of element which is not equals.
     *
     */
    template<typename LVector, typename RVector>
      typename std::enable_if<
      IsVector<LVector>::value && IsVector<RVector>::value,
      int
    >::type equalsReturnIndex (
      const LVector& lVector,
      const RVector& rVector,
      typename VectorTraits<LVector>::Scalar tolerance = NUMERICAL_ZERO_DIFFERENCE
      );


  } // namespace la
} // namespace tarch

#include "tarch/la/VectorVectorOperations.cpph"

#endif /* _TARCH_LA_VECTORVECTOROPERATIONS_H_ */
