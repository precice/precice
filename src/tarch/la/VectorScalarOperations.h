// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_LA_VECTORSCALAROPERATIONS_H_
#define _TARCH_LA_VECTORSCALAROPERATIONS_H_

#include "tarch/la/traits/IsVector.h"
#include "tarch/la/traits/VectorTraits.h"
#include "tarch/la/traits/DeduceScalar.h"
#include "tarch/utils/EnableIf.h"

namespace tarch {
  namespace la {

    /**
     * Multiplies every component of the vector with the scalar and assigns the
     * result to the vector.
     *
     * No temporary objects are created during the operation.
     */
    template<typename Vector>
      typename std::enable_if< IsVector<Vector>::value,
    Vector&>::type operator*= (
      Vector&                                      vector,
      const typename VectorTraits<Vector>::Scalar& scalar
    );

    /**
     * Divides every component of the vector by the scalar and assigns the
     * result to the vector.
     *
     * No temporary objects are created during the operation.
     */
    template<typename Vector>
      typename std::enable_if< IsVector<Vector>::value,
    Vector&>::type operator/= (
      Vector&                                      vector,
      const typename VectorTraits<Vector>::Scalar& scalar
    );

    /**
     * Adds every component of the vector to the scalar and assigns the
     * result to the vector.
     *
     * No temporary objects are created during the operation.
     */
    template<typename Vector>
      typename std::enable_if< IsVector<Vector>::value,
    Vector&>::type operator+= (
      Vector&                                      vector,
      const typename VectorTraits<Vector>::Scalar& scalar
    );

    /**
     * Subtracts the scalar from every component of the vector and assigns the
     * result to the vector.
     *
     * No temporary objects are created during the operation.
     */
    template<typename Vector>
      typename std::enable_if< IsVector<Vector>::value,
    Vector&>::type operator-= (
      Vector &                                      vector,
      const typename VectorTraits<Vector>::Scalar & scalar
    );

    /**
     * Multiplies every component of the vector with the scalar and returns the
     * result.
     *
     * A temporary vector is created during the operation and copied as result.
     */
    template<typename Vector>
      typename std::enable_if< IsVector<Vector>::value,
    Vector>::type operator* (
      const Vector&                                vector,
      const typename VectorTraits<Vector>::Scalar& scalar
    );

    /**
     * Divides every component of the vector by the scalar and returns the
     * result.
     *
     * A temporary vector is created during the operation and copied as result.
     */
    template<typename Vector>
      typename std::enable_if< IsVector<Vector>::value,
    Vector>::type operator/ (
      const Vector&                                vector,
      const typename VectorTraits<Vector>::Scalar& scalar
    );

    /**
     * Adds the scalar to every component of the vector and returns the
     * result.
     *
     * A temporary vector is created during the operation and copied as result.
     */
    template<typename Vector>
      typename std::enable_if< IsVector<Vector>::value,
    Vector>::type operator+ (
      const Vector&                                vector,
      const typename VectorTraits<Vector>::Scalar& scalar
    );

    /**
     * Subtracts the scalar from every component of the vector and returns the
     * result.
     *
     * A temporary vector is created during the operation and copied as result.
     */
    template<typename Vector>
      typename std::enable_if< IsVector<Vector>::value,
    Vector>::type operator- (
      const Vector&                                vector,
      const typename VectorTraits<Vector>::Scalar& scalar
    );

    /**
     * Multiplies every component of the vector with the scalar and returns the
     * result.
     *
     * A temporary vector is created during the operation and copied as result.
     */
    template<typename Vector>
      typename std::enable_if< IsVector<Vector>::value,
    Vector>::type operator* (
      const typename VectorTraits<Vector>::Scalar& scalar,
      const Vector&                                vector
    );

    /**
     * Adds the scalar to every component of the vector and returns the
     * result.
     *
     * A temporary vector is created during the operation and copied as result.
     */
    template<typename Vector>
      typename std::enable_if< IsVector<Vector>::value,
    Vector>::type operator+ (
      const typename VectorTraits<Vector>::Scalar& scalar,
      const Vector&                                vector
    );

    /**
     * Subtracts every component of the vector from the scalar and returns the
     * result.
     *
     * Attention: 1 - (v0, v1, v2) = (1 - v0, 1 - v1, 1 - v2).
     *
     * A temporary vector is created during the operation and copied as result.
     */
    template<typename Vector>
      typename std::enable_if< IsVector<Vector>::value,
    Vector>::type operator- (
      const typename VectorTraits<Vector>::Scalar& scalar,
      const Vector&                                vector
    );
    /**
     * Define operator %
     */
    template<typename Vector>
      typename std::enable_if< IsVector<Vector>::value,
    Vector>::type operator% (
        const Vector&                                vector,
        const typename VectorTraits<Vector>::Scalar& scalar
    );
  } // namespace la
} // namespace tarch

#include "tarch/la/VectorScalarOperations.cpph"

#endif /* _TARCH_LA_VECTORSCALAROPERATIONS_H_ */
