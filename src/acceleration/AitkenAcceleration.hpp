#pragma once

#include <Eigen/Core>
#include <algorithm>
#include <map>
#include <string>
#include <vector>
#include "acceleration/Acceleration.hpp"
#include "logging/Logger.hpp"

namespace precice {
namespace acceleration {

class AitkenAcceleration : public Acceleration {
public:
  AitkenAcceleration(
      double           initialRelaxationFactor,
      std::vector<int> dataIDs);

  virtual ~AitkenAcceleration() {}

  virtual std::vector<int> getDataIDs() const
  {
    return _dataIDs;
  }

  virtual void initialize(
      DataMap &cpldata);

  virtual void performAcceleration(
      DataMap &cpldata);

  virtual void iterationsConverged(
      DataMap &cpldata);

private:
  logging::Logger _log{"acceleration::AitkenAcceleration"};

  double _initialRelaxation;

  std::vector<int> _dataIDs;

  double _aitkenFactor;

  int _iterationCounter = 0;

  Eigen::VectorXd _residuals;
};
} // namespace acceleration
} // namespace precice
