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

  virtual ~ConstantRelaxationAcceleration() {}

  virtual std::vector<int> getDataIDs() const
  {
    return _dataIDs;
  }

  virtual void initialize(const DataMap &cplData);

  virtual void performAcceleration(const DataMap &cplData);

  virtual void iterationsConverged(const DataMap &cplData)
  {
  }

private:
  logging::Logger _log{"acceleration::ConstantRelaxationAcceleration"};

  double _relaxation;

  std::vector<int> _dataIDs;
};
} // namespace acceleration
} // namespace precice
