#pragma once

#include "tarch/la/Vector.h"
#include "tarch/la/DynamicVector.h"
#include "tarch/la/DynamicMatrix.h"
#include <Eigen/Dense>

namespace precice {
namespace utils {

typedef tarch::la::DynamicVector<double> DynVector;
typedef tarch::la::DynamicMatrix<double> DynMatrix;
typedef tarch::la::Vector<2,double>      Vector2D;
typedef tarch::la::Vector<3,double>      Vector3D;

std::string getTypeName(const DynVector& var);

std::string getTypeName(const Vector2D& var);

std::string getTypeName(const Vector3D& var);

const Eigen::VectorXd& delinearize (
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

/// @brief Provides mappings of indices for dimensions 2 and 3.
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

