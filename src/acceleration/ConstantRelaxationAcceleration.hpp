#pragma once

#include <algorithm>
#include <map>
#include <string>
#include <vector>
#include "acceleration/Acceleration.hpp"
#include "logging/Logger.hpp"

namespace precice {
namespace acceleration {

class ConstantRelaxationAcceleration : public Acceleration {
public:
  ConstantRelaxationAcceleration(
      double           relaxation,
      std::vector<int> dataIDs);

  virtual std::vector<int> getPrimaryDataIDs() const override final
  {
    return _dataIDs;
  }

  virtual void initialize(const DataMap &cplData) override;

  virtual void performAcceleration(DataMap &cplData) override;

  virtual void iterationsConverged(const DataMap &cplData) override
  {
    // function not needed in ConstantRelaxationAcceleration
  }

private:
  logging::Logger _log{"acceleration::ConstantRelaxationAcceleration"};

  double _relaxation;

  std::vector<int> _dataIDs;
};
} // namespace acceleration
} // namespace precice
