
#pragma once

#include <Eigen/Core>

namespace precice::time {

/// @brief Sample containing data values
struct Sample {
  Eigen::VectorXd values;
  Eigen::MatrixXd gradients;
};

} // namespace precice::time
