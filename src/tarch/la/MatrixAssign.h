// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_LA_MATRIXASSSIGN_H_
#define _TARCH_LA_MATRIXASSSIGN_H_

#include "tarch/la/traits/MatrixTraits.h"
#include "tarch/la/traits/IsMatrix.h"
#include "tarch/utils/EnableIf.h"
#include "tarch/Assertions.h"

namespace tarch {
  namespace la {
    template<typename Matrix> class MatrixAssign;

    /**
     * Returns a reinterpreted matrix to enable assignment of scalars and matrices.
     */
    template<typename Matrix>
      typename utils::EnableIf<IsMatrix<Matrix>::value,
      MatrixAssign<Matrix>&
    >::Type assign (
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
    typename utils::EnableIf< IsMatrix<RMatrix>::value,
    Matrix&
  >::Type operator= (const RMatrix& toAssign);
};

#include "tarch/la/MatrixAssign.cpph"

#endif /* _TARCH_LA_MATRIXASSSIGN_H_ */
