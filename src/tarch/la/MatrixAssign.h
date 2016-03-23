// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_LA_MATRIXASSSIGN_H_
#define _TARCH_LA_MATRIXASSSIGN_H_

#include "tarch/la/traits/MatrixTraits.h"
#include "tarch/la/traits/IsMatrix.h"
#include "utils/assertion.hpp"

namespace tarch {
  namespace la {
    template<typename Matrix> class MatrixAssign;

    /**
     * Returns a reinterpreted matrix to enable assignment of scalars and matrices.
     */
    template<typename Matrix>
      typename std::enable_if<IsMatrix<Matrix>::value,
      MatrixAssign<Matrix>&
    >::type assign (
      Matrix& matrix
    );
  }
}

/**
 * Reinterpreted matrix to enable assignment of scalars and matrices.
 */
template<typename Matrix>
class tarch::la::MatrixAssign
{
public:

  typedef MatrixTraits<Matrix> Traits;

  /**
   * Assigns a scalar to all components of the matrix.
   */
  Matrix& operator= (const typename Traits::Scalar& toAssign);

  /**
   * Assigns the components of a matrix to the matrix.
   */
  template<typename RMatrix>
    typename std::enable_if< IsMatrix<RMatrix>::value,
    Matrix&
  >::type operator= (const RMatrix& toAssign);
};

#include "tarch/la/MatrixAssign.cpph"

#endif /* _TARCH_LA_MATRIXASSSIGN_H_ */
