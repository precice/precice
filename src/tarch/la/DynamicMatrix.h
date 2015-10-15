// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_LA_DYNAMICMATRIX_H_
#define _TARCH_LA_DYNAMICMATRIX_H_

#include "tarch/la/MatrixAssign.h"
#include "tarch/la/MatrixAssignList.h"
#include "tarch/la/MatrixOperations.h"
#include "tarch/la/MatrixVectorOperations.h"
#include "tarch/la/MatrixMatrixOperations.h"
#include "Eigen/Core"
#include <vector>

namespace tarch {
  namespace la {
    template<typename Scalar> class DynamicMatrix;
  }
}

/**
 * Dynamic (i.e. runtime) sized matrix type, based on heap memory allocation.
 */
template<typename Scalar>
class tarch::la::DynamicMatrix
{
private:

  /**
   * Number of rows in the matrix.
   */
  int _rows;

  /**
   * Number of columns in the matrix.
   */
  int _cols;

  /**
   * Implementational basis for matrix element values.
   */
  std::vector<Scalar> _values;

public:

  /**
   * @brief Constructs an empty matrix.
   */
  DynamicMatrix ();

  /**
   * Constructs an uninitialized matrix of given size.
   */
  DynamicMatrix (int rows, int cols);

  /**
   * Constructs and initializes a matrix of given size.
   */
  DynamicMatrix (int rows, int cols, const Scalar& initialValue);

  /**
   * Returns the number of rows in the matrix.
   */
  int rows() const;

  /**
   * Returns the number of columns in the matrix.
   */
  int cols() const;

  /**
   * Returns the total number of elements in the matrix.
   */
  int size() const;

  /**
   * Returns matrix element at given row and column index.
   */
  Scalar& operator() (int rowIndex, int colIndex );

  /**
   * Returns const matrix element at given row and column index.
   */
  const Scalar& operator() (int rowIndex, int colIndex ) const;

  /// @brief Prints the matrix to stdout.
  void print() const;
  
  /// Converts to an Eigen::MatrixXd
  operator Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>() const;

  // No more methods here? They are all generic free methods now!
};

#include "tarch/la/DynamicMatrix.cpph"

#endif /* _TARCH_LA_DYNAMICMATRIX_H_ */
