#ifndef _TARCH_LA_MATRIXMATRIXOPERATIONS_H_
#define _TARCH_LA_MATRIXMATRIXOPERATIONS_H_

#include "tarch/la/traits/IsMatrix.h"
#include "tarch/la/traits/MatrixTraits.h"
#include "tarch/la/Scalar.h"

namespace tarch {
  namespace la {

    /**
     * Performs a matrix x matrix multiplication.
     *
     * The resulting matrix has to be created outside and given as a parameter.
     */
    template<typename LMatrix, typename RMatrix, typename ResultMatrix>
      typename std::enable_if<
      IsMatrix<LMatrix>::value && IsMatrix<RMatrix>::value && IsMatrix<ResultMatrix>::value,
      void
    >::type multiply (
      const LMatrix& lMatrix,
      const RMatrix& rMatrix,
      ResultMatrix&  result);

    /**
     * Bitwise comparison of the components of two matrices on equality.
     */
    template<typename LMatrix, typename RMatrix>
      typename std::enable_if<
      IsMatrix<LMatrix>::value && IsMatrix<RMatrix>::value,
      bool
    >::type operator== (
      const LMatrix& lMatrix,
      const RMatrix& rMatrix);

    /**
     * Compares to matrices on equality by means of a numerical accuracy.
     */
    template<typename LMatrix, typename RMatrix>
      typename std::enable_if<
      IsMatrix<LMatrix>::value && IsMatrix<RMatrix>::value /*&& EqualScalars<LMatrix,RMatrix>::value*/,
      bool
    >::type equals (
      const LMatrix&                         lMatrix,
      const RMatrix&                         rMatrix,
      typename MatrixTraits<LMatrix>::Scalar tolerance = NUMERICAL_ZERO_DIFFERENCE
    );
    /**
     * Adds rVector to lVector.
     *
     * A temporary vector is created and copied to store return back the result.
     */
    template<typename LMatrix, typename RMatrix>
      typename std::enable_if<
      IsMatrix<LMatrix>::value && IsMatrix<RMatrix>::value,
      LMatrix
    >::type operator+ (
      const LMatrix& lMatrix,
      const RMatrix& rMatrix
      );
    /**
     * Return Index of element which is not equals.
     *
     */
    template<typename LMatrix, typename RMatrix>
      typename std::enable_if<
      IsMatrix<LMatrix>::value && IsMatrix<RMatrix>::value,
      int
    >::type equalsReturnIndex (
      const LMatrix& lMatrix,
      const RMatrix& rMatrix,
      typename MatrixTraits<LMatrix>::Scalar tolerance = NUMERICAL_ZERO_DIFFERENCE
      );

  }
} // namespace tarch, la

#include "tarch/la/MatrixMatrixOperations.cpph"

#endif /* _TARCH_LA_MATRIXMATRIXOPERATIONS_H_ */
