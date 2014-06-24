// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_LA_GRAMSCHMIDT_H_
#define _TARCH_LA_GRAMSCHMIDT_H_

#include "tarch/la/traits/IsMatrix.h"
#include "tarch/la/traits/MatrixTraits.h"
#include "tarch/la/traits/IsVector.h"
#include "tarch/la/traits/VectorTraits.h"
#include "tarch/la/VectorOperations.h"
#include "tarch/la/VectorVectorOperations.h"
#include "tarch/la/VectorScalarOperations.h"
#include "tarch/la/VectorAssign.h"
#include "tarch/utils/EnableIf.h"
#include "tarch/Assertions.h"

namespace tarch {
  namespace la {

  /**
   * Produces an economy-size QR decomposition of a matrix A, A is changed.
   *
   * Preconditions:
   * - The column vectors of A have to be linear independent. Otherwise, see
   *   postconditions "R".
   * - If A has n columns and m rows, then:
   *   * n <= m for A
   *   * Q has to have same size than A
   *   * R has to be a n * n matrix
   *
   * Postconditions:
   * - A is changed
   * - Q holds an orthonormal basis of range(A) in its column vectors.
   * - R is an upper triangular matrix, with all main diagonal element unequal
   *   to zero iff A has full rank.
   *
   * The algorithm is taken from the book Numerical Mathematics, Second Edition,
   * by Alfio Querteroni, Riccardo Sacco, and Fausto Saleri, Springer, 2007,
   * pages 85-87.
   */
  template<typename Matrix>
    typename utils::EnableIf<IsMatrix<Matrix>::value,
    void
  >::Type modifiedGramSchmidt (
    Matrix& A,
    Matrix& Q,
    Matrix& R
  );

  } // namespace la
} // namespace tarch

#include "tarch/la/GramSchmidt.cpph"

#endif /* PRECICE_LA_GRAMSCHMIDT_HPP */
