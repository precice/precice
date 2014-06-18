// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_LA_TRAITS_ISMATRIX_H_
#define _TARCH_LA_TRAITS_ISMATRIX_H_

namespace tarch {
   namespace la {
      template<int Rows, int Cols, typename Scalar> class Matrix;
      template<typename Scalar> class DynamicMatrix;
      template<typename Scalar> class DynamicColumnMatrix;
      template<typename Matrix> class TransposedMatrix;
   }
}


namespace tarch {
namespace la {

template<typename Matrix>
struct IsMatrix
{
   static const int value = false;
};

template<int Rows, int Cols, typename Scalar>
struct IsMatrix<Matrix<Rows,Cols,Scalar> >
{
   static const int value = true;
};

template<typename Scalar>
struct IsMatrix<DynamicMatrix<Scalar> >
{
   static const int value = true;
};

template<typename Scalar>
struct IsMatrix<DynamicColumnMatrix<Scalar> >
{
   static const int value = true;
};

template<typename Matrix>
struct IsMatrix<TransposedMatrix<Matrix> >
{
   static const int value = true;
};

}} // namespace tarch, la

#endif /* _TARCH_LA_TRAITS_ISMATRIX_H_ */
