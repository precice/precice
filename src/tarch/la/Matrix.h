// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_LA_MATRIX_H_
#define _TARCH_LA_MATRIX_H_

#include "tarch/la/MatrixAssign.h"
#include "tarch/la/MatrixAssignList.h"
#include "tarch/la/MatrixOperations.h"
#include "tarch/la/MatrixVectorOperations.h"
#include "tarch/la/MatrixMatrixOperations.h"
#include "tarch/la/MatrixScalarOperations.h"

namespace tarch {
  namespace la {
    template<int Rows, int Cols, typename Scalar> class Matrix;
  }
}

/**
 * Static (i.e. compiletime) sized matrix type.
 */
template<int Rows, int Cols, typename Scalar>
class tarch::la::Matrix
{
private:

  /**
   * Values of the matrix components.
   */
  Scalar _values[Rows*Cols];

public:

  /**
   * Constructs a non-initialized matrix.
   */
  Matrix ();

  /**
   * Constructs an initialized matrix.
   */
  Matrix (const Scalar& initialValue);

  /**
   * Returns the number of rows in the matrix.
   */
  int rows() const;

  /**
   * Returns the number of columns in the matrix.
   */
  int cols() const;

  /**
   * Returns the number of total elements in the matrix.
   */
  int size() const;

  /**
   * Returns element at given row and column index (from 0..size-1).
   */
  Scalar & operator() (
    int rowIndex,
    int colIndex );

  /**
   * Returns const element at given row and column index (from 0..size-1).
   */
  const Scalar & operator() (
    int rowIndex,
    int colIndex ) const;


  // No more methods here? They are all generic free methods now!
};

#include "tarch/la/Matrix.cpph"

#endif /* _TARCH_LA_MATRIX_H_ */
