
#pragma once

#include <Eigen/Core>

namespace precice::time {

/// @brief Sample of a \ref mesh::Data on a \ref mesh::Mesh
struct Sample {
  /// Constructs an invalid empty Sample
  Sample() = default;

  /// Constructs an empty Sample of a given data dimensionality
  explicit Sample(int dims) noexcept
      : dataDims(dims) {}

  /// Constructs a Sample of given data dimensionality and size without gradients
  Sample(int dims, int dataCount)
      : dataDims(dims), values(dims * dataCount) {}

  /// Constructs a Sample of given data dimensionality and data values
  Sample(int dims, Eigen::VectorXd inValues)
      : dataDims(dims), values(inValues) {}

  /// Constructs a Sample of given data dimensionality, data values, and data gradients
  Sample(int dims, Eigen::VectorXd inValues, Eigen::MatrixXd inGradients)
      : dataDims(dims), values(inValues), gradients(inGradients) {}

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
  int dataDims{-1};

  /// The data values linearised
  /// @todo Change to matrix so that values.col(i) gets the value at vertex i
  Eigen::VectorXd values;

  /// The gradients of the data. Use gradients.col(i) to get the gradient at vertex i
  Eigen::MatrixXd gradients;
};

} // namespace precice::time
