// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_LA_VECTOR_COMPARE_H_
#define _TARCH_LA_VECTOR_COMPARE_H_


#include "tarch/la/Vector.h"


namespace tarch {
  namespace la {
    template<int N>
    struct VectorCompare;
  }
}

/**
 * Comparison operator for boolean vectors.
 *
 * If you wanna store sets of vectors (e.g. in a map), the STL libraries are by
 * default falling back to a bit-wise comparison. This might not be of use if
 * you manage sets/maps of double vectors. In this case, you wanna the container
 * take into account machine precision. In such a case, declare the map e.g. as
 * follows:
 *  \code
  std::map<tarch::la::Vector<DIMENSIONS,double> , int, tarch::la::VectorCompare<DIMENSIONS>  >  _vertex2IndexMap;
 \endcode
 * @author Tobias Weinzierl
 */
template<int N>
struct tarch::la::VectorCompare {
  bool operator()(
    const Vector<N,double>& left,
    const Vector<N,double>& right
  ) const {
    return firstGreater(right, left);
  }
};


#endif
