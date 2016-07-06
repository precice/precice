#ifndef _TARCH_LA_TRANSPOSEDMATRIX_H_
#define _TARCH_LA_TRANSPOSEDMATRIX_H_

#include "tarch/la/MatrixAssign.h"
#include "tarch/la/MatrixAssignList.h"
#include "tarch/la/MatrixOperations.h"
#include "tarch/la/MatrixVectorOperations.h"
#include "tarch/la/MatrixMatrixOperations.h"
#include "tarch/la/traits/IsMatrix.h"
#include "tarch/la/traits/MatrixTraits.h"

namespace tarch {
  namespace la {
    template<typename Matrix> class TransposedMatrix;

    /**
     * Returns a reference to the given matrix with transposed traits semantics.
     * This operation has zero performance or memory overhead if properly inlined.
     */
    template<typename Matrix>
      typename std::enable_if<IsMatrix<Matrix>::value,
      TransposedMatrix<Matrix>&
    >::type transpose (Matrix& matrix);
  }
}

template<typename Matrix>
class tarch::la::TransposedMatrix
{
public:

  typedef typename MatrixTraits<Matrix>::Scalar Scalar;

  int rows() const;

  int cols() const;

  int size() const;

  Scalar & operator() (
    int rowIndex,
    int colIndex );

  const Scalar & operator() (
    int rowIndex,
    int colIndex ) const;
};

#include "tarch/la/TransposedMatrix.cpph"

#endif /* _TARCH_LA_TRANSPOSEDMATRIX_H_ */
