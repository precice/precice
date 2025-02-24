#pragma once

#include <algorithm>
#include <map>
#include <string>
#include <vector>
#include "acceleration/Acceleration.hpp"
#include "logging/Logger.hpp"

namespace precice::acceleration {

class ConstantRelaxationAcceleration : public Acceleration {
public:
  ConstantRelaxationAcceleration(
      double           relaxation,
      std::vector<int> dataIDs);

  std::vector<int> getPrimaryDataIDs() const final
  {
    return _dataIDs;
  }

  void initialize(const DataMap &cplData) override;

  void performAcceleration(DataMap &cplData, double windowStart, double windowEnd) override;

  void iterationsConverged(const DataMap &cplData, double windowStart) override
  {
    // function not needed in ConstantRelaxationAcceleration
  }

private:
  logging::Logger _log{"acceleration::ConstantRelaxationAcceleration"};

  double _relaxation;

  std::vector<int> _dataIDs;
};
} // namespace precice::acceleration
