// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_LA_MATRIXSCALAROPERATIONS_H_
#define _TARCH_LA_MATRIXSCALAROPERATIONS_H_

#include "tarch/la/traits/IsMatrix.h"
#include "tarch/la/traits/MatrixTraits.h"
#include "tarch/la/traits/DeduceScalar.h"

namespace tarch {
  namespace la {
    /**
     * Multiplies every component of the vector with the scalar and assigns the
     *  result to the vector.
     *
     * No temporary objects are created during the operation.
     */
     template<typename Matrix>
       typename std::enable_if< IsMatrix<Matrix>::value,
     Matrix&>::type operator*= (
       Matrix&                                      matrix,
       const typename MatrixTraits<Matrix>::Scalar& scalar
    );
     /**
      * Multiplies every component of the vector with the scalar and assigns the
      *  result to the vector.
      *
      * No temporary objects are created during the operation.
      */
      template<typename Matrix>
        typename std::enable_if< IsMatrix<Matrix>::value,
      Matrix>::type operator* (
        const Matrix&                                matrix,
        const typename MatrixTraits<Matrix>::Scalar& scalar
     );
  }
}

#include "tarch/la/MatrixScalarOperations.cpph"

#endif
