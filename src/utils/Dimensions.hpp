#pragma once

#include <Eigen/Core>

namespace precice {
namespace utils {

const Eigen::VectorXd &delinearize(
    int toDelinearize,
    int dimensions);

template <typename VECTOR_T>
int linearize(
    const VECTOR_T &toLinearize)
{
  int index = 0;
  for (int dim = 0; dim < toLinearize.size(); dim++) {
    if (toLinearize(dim) > 0.0) {
      index += (int) std::pow(2.0, dim);
    }
  }
  return index;
}

/// Provides mappings of indices for dimensions 2 and 3.
template <int dimension>
struct IndexMaps {
};

template <>
struct IndexMaps<2> {
  static const int CUBOID_EDGE_VERTICES[4][2];
};

template <>
struct IndexMaps<3> {
  static const int CUBOID_FACE_VERTICES[6][4];
  static const int CUBOID_FACE_EDGES[6][4];
  static const int CUBOID_EDGE_VERTICES[12][2];
};

} // namespace utils
} // namespace precice
