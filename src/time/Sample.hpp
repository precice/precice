
#pragma once

#include <Eigen/Core>

namespace precice::time {

/// @brief Sample containing data values
struct Sample {
  Sample() = default;
  explicit Sample(int dims) noexcept
      : dataDims(dims) {}
  Sample(int dims, int dataCount)
      : dataDims(dims), values(dims * dataCount) {}
  Sample(int dims, Eigen::VectorXd inValues)
      : dataDims(dims), values(inValues) {}
  Sample(int dims, Eigen::VectorXd inValues, Eigen::MatrixXd inGradients)
      : dataDims(dims), values(inValues), gradients(inGradients) {}

  Sample(const Sample &) = default;
  Sample(Sample &&)      = default;

  Sample &setZero()
  {
    values.setZero();
    gradients.setZero();
    return *this;
  }

  Sample &operator=(const Sample &) = default;
  Sample &operator=(Sample &&) = default;

  int             dataDims{-1};
  Eigen::VectorXd values;
  Eigen::MatrixXd gradients;
};

} // namespace precice::time
