// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_UTILS_DIMENSIONS_HPP_
#define PRECICE_UTILS_DIMENSIONS_HPP_

#include "tarch/la/Vector.h"
#include "tarch/la/DynamicVector.h"

namespace precice {
namespace utils {

/**
 * @brief Compile-time computation of base "to the power of" exp
 *
 * Exp and base must be of type int
 */
template<int base, int exp> struct CompilePower {
   enum { VALUE = base * CompilePower<base, exp-1>::VALUE };
};
template< int base > struct CompilePower<base, 1> {
   enum { VALUE = base };
};

typedef tarch::la::DynamicVector<double> DynVector;
typedef tarch::la::Vector<2,double>      Vector2D;
typedef tarch::la::Vector<3,double>      Vector3D;

std::string getTypeName(const DynVector& var);

std::string getTypeName(const Vector2D& var);

std::string getTypeName(const Vector3D& var);

const utils::DynVector& delinearize (
  int toDelinearize,
  int dimensions );

template<typename VECTOR_T>
int linearize
(
  const VECTOR_T& toLinearize )
{
  int index = 0;
  for (int dim=0; dim < toLinearize.size(); dim++) {
    if (toLinearize(dim) > 0.0) {
      index += (int) std::pow(2.0, dim);
    }
  }
  return index;
}

/**
 * @brief Provides mappings of indices for dimensions 2 and 3.
 */
template< int dimension > struct IndexMaps {};

template<> struct IndexMaps<2>
{
   static const int CUBOID_EDGE_VERTICES[4][2];
};

template<> struct IndexMaps<3>
{
   static const int CUBOID_FACE_VERTICES[6][4];
   static const int CUBOID_FACE_EDGES[6][4];
   static const int CUBOID_EDGE_VERTICES[12][2];
};

}} // namespace precice, utils

#endif /* PRECICE_UTILS_DIMENSIONS_HPP_ */
