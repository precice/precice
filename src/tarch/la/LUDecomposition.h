// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_LA_LUDECOMPOSITION_H_
#define _TARCH_LA_LUDECOMPOSITION_H_

#include "tarch/la/traits/IsMatrix.h"
#include "tarch/la/traits/MatrixTraits.h"
#include "tarch/la/traits/IsVector.h"
#include "tarch/la/traits/VectorTraits.h"
#include "tarch/utils/EnableIf.h"
#include "tarch/Assertions.h"

namespace tarch {
  namespace la {

//    /**
//     * Performs an inplace LU-decomposition of the square matrix A.
//     */
//    template<typename Matrix>
//      typename utils::EnableIf<IsMatrix<Matrix>::value,
//      Matrix&
//    >::Type luCompact (
//      Matrix& A
//    ) {
//      DynamicVector<double> pivots;
//      return lu(A,pivots);
//    }

    /**
     * Performs an inplace LU-decomposition of the square matrix A. Adds pivot
     * line indices.
     */
    template<typename Matrix, typename Vector>
      typename utils::EnableIf<IsMatrix<Matrix>::value && IsVector<Vector>::value,
      Matrix&
    >::Type lu (
      Matrix& A,
      Vector& pivots
    ) {
      typedef MatrixTraits<Matrix> MT;
      typedef VectorTraits<Vector> VT;
      assertion2 (MT::rows(A) == MT::cols(A), MT::rows(A), MT::cols(A));
      int n = MT::rows(A); // quadratic matrix size
      assertion2 (n == VT::size(pivots), n, VT::size(pivots));

      // Init pivots
      //for (int i=0; i < n; i++){
      //  VT::elem(i,pivots) = i;
      //}

      // Perform LU-decomposition
      for (int k=0; k < n; k++){
        //std::cout << "Row " << k << ": " << A << std::endl;
        // Determine line with max pivot
        int maxIndex = k;
        typename MT::Scalar max = std::abs(MT::celem(k,k,A));
        for (int i=k+1; i < n; i++){
          typename MT::Scalar current = std::abs(MT::celem(i,k,A));
          if (current > max){
            maxIndex = i;
            max = current;
          }
        }
        //assertion1 (greater(max,0.0), max);
        //std::cout << "Exchange row " << k << " with row " << maxIndex << std::endl;

        // Change pivots
        //int temp = VT::celem(k,pivots);
        //VT::elem(k,pivots) = VT::celem(maxIndex,pivots);
        //VT::elem(maxIndex,pivots) = temp;
        VT::elem(k,pivots) = maxIndex;

        // Exchange lines
        for (int j=0; j < n; j++){
          typename MT::Scalar temp = MT::celem(k,j,A);
          MT::elem(k,j,A) = MT::elem(maxIndex,j,A);
          MT::elem(maxIndex,j,A) = temp;
        }

        // Compute scaling elements
        for (int i=k+1; i < n; i++){
          MT::elem(i,k,A) /= MT::celem(k,k,A);
        }

        // Subtract contributions from each line
        for (int i=k+1; i < n; i++){
          for (int j=k+1; j < n; j++){
            MT::elem(i,j,A) -= MT::celem(i,k,A) * MT::celem(k,j,A);
          }
        }
      }
      return A;
    }

  } // namespace la
} // namespace tarch

#endif /* _TARCH_LA_LUDECOMPOSITION_H_ */
