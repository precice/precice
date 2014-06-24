// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_LA_MATRIXVECTOROPERATIONS_H_
#define _TARCH_LA_MATRIXVECTOROPERATIONS_H_

#include "tarch/la/traits/IsMatrix.h"
#include "tarch/la/traits/IsVector.h"
#include "tarch/la/traits/MatrixTraits.h"
#include "tarch/la/traits/VectorTraits.h"
#include "tarch/utils/EnableIf.h"

namespace tarch {
  namespace la {

    /**
     * Performs a matrix x vector multiplication.
     *
     * The result vector has to be created outside and given as a parameter.
     */
    template<typename Matrix, typename Vector, typename Result>
      typename utils::EnableIf<
      IsMatrix<Matrix>::value && IsVector<Vector>::value && IsVector<Result>::value,
      Result&
    >::Type multiply (
      const Matrix& matrix,
      const Vector& vector,
      Result&       result );

    /**
     * Performs a matrix x vector multiplication.
     *
     * The result vector is created inside, multiply is used to execute the
     * multiplication.
     */
    template<typename Matrix, typename Vector>
      typename utils::EnableIf<
      IsMatrix<Matrix>::value && IsVector<Vector>::value,
      Vector
    >::Type operator* (
      const Matrix& matrix,
      const Vector& vector );

    /**
     * Solvers the 3 by 3 linear system: matrix * result = rhs.
     */
    template<typename Matrix, typename Vector, typename Result>
      typename utils::EnableIf<
      IsMatrix<Matrix>::value && IsVector<Vector>::value && IsVector<Vector>::value,
      Result&
    >::Type solveSystem3x3 (
      const Matrix& matrix,
      const Vector& rhs,
      Result&       result );

    /**
     * Performs a forward-substitution of the system: L * x = rhs.
     *
     * @param matrix A lower trianguler matrix, diagonal values are assumed to be
     *        1.0. Upper triangular values and diagonal can be arbitrary.
     * @param rhs    The righthandside of the equation system.
     * @param x      Unknown vector to be solved for.
     */
    template<typename Matrix, typename Vector>
      typename utils::EnableIf< IsMatrix<Matrix>::value && IsVector<Vector>::value,
      void
    >::Type forwardSubstitution (
      const Matrix& matrix,
      const Vector& rhs,
      Vector&       x );

    /**
     * Performs a back-substitution of the system: U * x = rhs.
     *
     * @param matrix An upper trianguler matrix (lower triangle can be arbitrary).
     * @param rhs    The righthandside of the equation system.
     * @param x      Unknown vector to be solved for.
     */
    template<typename Matrix, typename Vector>
      typename utils::EnableIf< IsMatrix<Matrix>::value && IsVector<Vector>::value,
      void
    >::Type backSubstitution (
      const Matrix& matrix,
      const Vector& rhs,
      Vector&       x );

  }
}

#include "MatrixVectorOperations.cpph"

#endif /* _TARCH_LA_MATRIXVECTOROPERATIONS_H_ */
