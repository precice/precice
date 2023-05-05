
#pragma once

#include <Eigen/Core>

namespace precice::time {

/// @brief Sample containing data values
struct Sample {
  Eigen::VectorXd values;
  Eigen::MatrixXd gradients; // @todo rename to gradientValues? See Data and CouplingData
};

} // namespace precice::time
