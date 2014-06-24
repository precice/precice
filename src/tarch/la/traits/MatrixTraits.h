// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_LA_TRAITS_MATRIXTRAITS_H_
#define _TARCH_LA_TRAITS_MATRIXTRAITS_H_

namespace tarch {
  namespace la {
    template<int Rows, int Cols, typename SCALAR> class Matrix;
    template<typename SCALAR> class DynamicMatrix;
    template<typename Scalar> class DynamicColumnMatrix;
    template<typename Matrix> class TransposedMatrix;
  }
}

namespace tarch {
namespace la {

template<typename Matrix>
struct MatrixTraits {};

template<int Rows, int Cols, typename SCALAR>
struct MatrixTraits<Matrix<Rows,Cols,SCALAR> >
{
  typedef Matrix<Rows,Cols,SCALAR> ThisMatrix;
  typedef SCALAR Scalar;

  static int rows (const ThisMatrix& matrix)
  {
    return Rows;
  }

  static int cols (const ThisMatrix& matrix)
  {
    return Cols;
  }

  static SCALAR& elem (
    int      rowIndex,
    int      colIndex,
    ThisMatrix& matrix
  ) {
    return matrix(rowIndex, colIndex);
  }

  static const SCALAR& celem (
    int            rowIndex,
    int            colIndex,
    const ThisMatrix& matrix
  ) {
    return matrix(rowIndex, colIndex);
  }
};

template<typename SCALAR>
struct MatrixTraits<DynamicMatrix<SCALAR> >
{
  typedef DynamicMatrix<SCALAR> ThisMatrix;
  typedef SCALAR Scalar;

  static int rows (const ThisMatrix& matrix) {
    return matrix.rows();
  }

  static int cols (const ThisMatrix& matrix) {
    return matrix.cols();
  }

  static SCALAR& elem (
    int      rowIndex,
    int      colIndex,
    ThisMatrix& matrix
  ) {
    return matrix(rowIndex, colIndex);
  }

  static const SCALAR& celem (
    int            rowIndex,
    int            colIndex,
    const ThisMatrix& matrix )
  {
    return matrix(rowIndex, colIndex);
  }
};

template<typename SCALAR>
struct MatrixTraits<DynamicColumnMatrix<SCALAR> >
{
  typedef DynamicColumnMatrix<SCALAR> ThisMatrix;
  typedef SCALAR Scalar;

  static int rows (const ThisMatrix& matrix) {
    return matrix.rows();
  }

  static int cols (const ThisMatrix& matrix) {
    return matrix.cols();
  }

  static SCALAR& elem (
    int      rowIndex,
    int      colIndex,
    ThisMatrix& matrix
  ) {
    return matrix(rowIndex, colIndex);
  }

  static const SCALAR& celem (
    int            rowIndex,
    int            colIndex,
    const ThisMatrix& matrix )
  {
    return matrix(rowIndex, colIndex);
  }
};

template<typename Matrix>
struct MatrixTraits<TransposedMatrix<Matrix> >
{
  typedef TransposedMatrix<Matrix> ThisMatrix;
  typedef typename TransposedMatrix<Matrix>::Scalar Scalar;

  static int rows (const ThisMatrix& matrix) {
    return matrix.rows();
  }

  static int cols (const ThisMatrix& matrix)
  {
    return matrix.cols();
  }

  static Scalar& elem (
    int      rowIndex,
    int      colIndex,
    ThisMatrix& matrix
  ) {
    return matrix(rowIndex, colIndex);
  }

  static const Scalar& celem (
    int            rowIndex,
    int            colIndex,
    const ThisMatrix& matrix )
  {
    return matrix(rowIndex, colIndex);
  }
};

}} // namespace tarch, la

#endif /* _TARCH_LA_TRAITS_MATRIXTRAITS_H_ */
