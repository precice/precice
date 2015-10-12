// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_LA_MATRIXOPERATIONS_H_
#define _TARCH_LA_MATRIXOPERATIONS_H_

#include "tarch/la/traits/IsMatrix.h"
#include "tarch/la/traits/MatrixTraits.h"
#include "tarch/utils/EnableIf.h"
#include <sstream>
#include <cmath>

namespace tarch {
  namespace la {

    /**
     * Computes the determinant of a 3 by 3 matrix.
     */
    template<typename Matrix>
      typename utils::LazyEnableIf< IsMatrix<Matrix>::value,
      utils::LazyType<typename MatrixTraits<Matrix>::Scalar>
    >::Type det3x3 ( const Matrix& matrix );


    /**
     * Computes the sum of all entries of the matrix.
     */
    template<typename Matrix>
      typename utils::LazyEnableIf< IsMatrix<Matrix>::value,
      utils::LazyType<typename MatrixTraits<Matrix>::Scalar>
    >::Type sum (const Matrix& matrix);
    
    /**
     * Computes the frobenius norm of the matrix
     */
    template<typename Matrix>
      typename utils::LazyEnableIf< IsMatrix<Matrix>::value,
      utils::LazyType<typename MatrixTraits<Matrix>::Scalar>
    >::Type frobeniusNorm (const Matrix& matrix);



    /**
     * Computes the square root of every component of the matrix.
     */
    template<typename Matrix>
      typename utils::EnableIf< IsMatrix<Matrix>::value,
      Matrix>::Type sqrt (
      const Matrix&                               matrix);
  }
}

template<typename Matrix>
  typename tarch::utils::EnableIf< tarch::la::IsMatrix<Matrix>::value,
  std::ostream&
>::Type operator<< (std::ostream& os, const Matrix& matrix);

#include "tarch/la/MatrixOperations.cpph"

#endif /* _TARCH_LA_MATRIXOPERATIONS_H_ */
