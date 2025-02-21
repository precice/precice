
#pragma once

#include <Eigen/Core>

#include "utils/assertion.hpp"

namespace precice::time {

/** Sample of a \ref mesh::Data on a \ref mesh::Mesh
 *
 * A \ref Sample encapsulates user-provided data of a mesh including values and gradients.
 * It is also aware of the dimensionality of the data, which is important for empty meshes.
 */
struct Sample {
  /// Constructs an empty Sample of a given data dimensionality
  explicit Sample(int dims) noexcept
      : dataDims(dims) {}

  /// Constructs a Sample of given data dimensionality and size without gradients
  Sample(int dims, int dataCount)
      : dataDims(dims), values(dims * dataCount) {}

  /// Constructs a Sample of given data and mesh dimensionality, and size with gradients
  Sample(int dataDims, int nVertices, int meshDims)
      : dataDims(dataDims), values(nVertices * dataDims), gradients(meshDims, nVertices * dataDims)
  {
  }

  /// Constructs a Sample of given data dimensionality and data values
  Sample(int dims, Eigen::VectorXd inValues)
      : dataDims(dims), values(std::move(inValues))
  {
    PRECICE_ASSERT(dataDims > 0);
  }

  /// Constructs a Sample of given data dimensionality, data values, and data gradients
  Sample(int dims, Eigen::VectorXd inValues, Eigen::MatrixXd inGradients)
      : dataDims(dims), values(std::move(inValues)), gradients(std::move(inGradients))
  {
    PRECICE_ASSERT(dataDims > 0);
    PRECICE_ASSERT(gradients.size() == 0 || gradients.cols() == values.size());
  }

  Sample(const Sample &) = default;
  Sample(Sample &&)      = default;

  Sample &operator=(const Sample &) = default;
  Sample &operator=(Sample &&) = default;

  /// Sets values and gradients to zero
  Sample &setZero()
  {
    values.setZero();
    gradients.setZero();
    return *this;
  }

  /// The dimensionality of the data
  int dataDims;

  /// The data values linearised
  /// @todo Change to matrix so that values.col(i) gets the value at vertex i
  Eigen::VectorXd values;

  /// The gradients of the data. Use gradients.col(i) to get the gradient at vertex i
  Eigen::MatrixXd gradients;
};

} // namespace precice::time
