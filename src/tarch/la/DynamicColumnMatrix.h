#ifndef _TARCH_LA_DYNAMICCOLUMNMATRIX_H_
#define _TARCH_LA_DYNAMICCOLUMNMATRIX_H_

#include "tarch/la/MatrixAssign.h"
#include "tarch/la/MatrixAssignList.h"
#include "tarch/la/MatrixOperations.h"
#include "tarch/la/MatrixVectorOperations.h"
#include "tarch/la/MatrixMatrixOperations.h"
#include "tarch/la/MatrixScalarOperations.h"
#include "tarch/la/DynamicVector.h"
#include <deque>

namespace tarch {
  namespace la {
    template<typename Scalar> class DynamicColumnMatrix;
  }
}

/**
 * Dynamic (i.e. runtime) sized matrix type, based on heap memory allocation.
 */
template<typename Scalar>
class tarch::la::DynamicColumnMatrix
{
private:

  /**
   * Number of rows in the matrix.
   */
  int _rows;

  /**
   * Matrix column vectors.
   */
  std::deque<DynamicVector<Scalar> > _columnVectors;

public:

  /**
   * Constructs an empty matrix.
   */
  DynamicColumnMatrix();

  /**
   * Constructs an uninitialized matrix of given size.
   */
  DynamicColumnMatrix(int rows, int cols);

  /**
   * Constructs and initializes a matrix of given size.
   */
  DynamicColumnMatrix(int rows, int cols, const Scalar& initialValue);

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
   * Returns a column of the matrix. Efficient, due to column based layout.
   */
  const DynamicVector<Scalar>& column(int colIndex) const;

  /**
   * Returns a column of the matrix. Efficient, due to column based layout.
   */
  DynamicVector<Scalar>& column(int colIndex);

  /**
   * Appends a column to front of matrix.
   *
   * Efficient due to column based storage layout.
   */
  void appendFront(const DynamicVector<Scalar>& columnVector);

  /**
   * Appends a column to back of matrix.
   *
   * Efficient due to column based storage layout.
   */
  void append(const DynamicVector<Scalar>& columnVector);

  /**
   * Appends a matrix column-wise to back of matrix.
   *
   * Efficient due to column based storage layout.
   */
  void append(const DynamicColumnMatrix<Scalar>& matrix);

  /**
   * Shifts all entries one to the right, last column entries are lost. Sets
   * entries of given columnVector at first column of matrix.
   *
   * Efficient due to column based storage layout.
   */
  void shiftSetFirst(const DynamicVector<Scalar>& columnVector);

  /**
   * Removes a column from the matrix. This operation involves reordering of
   * a std::Vector, which is not very efficient.
   */
  void remove(int colIndex);

  /**
   * Removes all entries from the matrix.
   */
  void clear();

  /**
   * Returns matrix element at given row and column index.
   */
  Scalar& operator()(int rowIndex, int colIndex);

  /**
   * Returns const matrix element at given row and column index.
   */
  const Scalar& operator()(int rowIndex, int colIndex) const;

  void printm(const char* filename) const;

  // No more methods here? They are all generic free methods now!
};

#include "tarch/la/DynamicColumnMatrix.cpph"

#endif /* _TARCH_LA_DYNAMICCOLUMNMATRIX_H_ */
