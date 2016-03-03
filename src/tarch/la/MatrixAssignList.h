// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_LA_MATRIXASSSIGNLIST_H_
#define _TARCH_LA_MATRIXASSSIGNLIST_H_

#include "tarch/la/traits/MatrixTraits.h"
#include "tarch/la/traits/IsMatrix.h"
#include "tarch/Assertions.h"

namespace tarch {
  namespace la {
    template<typename Matrix> class MatrixAssignList;

    /**
     * Returns a wrapper around a matrix to enable comma separated list assignment.
     */
    template<typename Matrix>
      typename std::enable_if<IsMatrix<Matrix>::value,
      MatrixAssignList<Matrix>
    >::type assignList ( Matrix& matrix );
  }
}

/**
 * Wrapper around a matrix to enable comma separated list assignment.
 */
template<typename Matrix>
class tarch::la::MatrixAssignList
{
private:

  Matrix& _matrix;
  int _rowIndex;
  int _colIndex;

public:

  typedef MatrixTraits<Matrix> Traits;

  /**
   * Constructor.
   */
  MatrixAssignList (Matrix & matrix);

  /**
   * Destructor, asserts that the right amount of values has been assigned.
   */
  ~MatrixAssignList();

  /**
   * Assigns the first value.
   */
  MatrixAssignList<Matrix> & operator= (const typename Traits::Scalar & toAssign);

  /**
   * Assigns all other values after operator=() has been invoked.
   */
  MatrixAssignList<Matrix> & operator, (const typename Traits::Scalar & toAssign);
};

#include "tarch/la/MatrixAssignList.cpph"

#endif /* _TARCH_LA_MATRIXASSSIGNLIST_H_ */
