#pragma once

#include <Eigen/Core>
#include <algorithm>
#include <map>
#include <string>
#include <vector>
#include "acceleration/Acceleration.hpp"
#include "acceleration/impl/Preconditioner.hpp"
#include "acceleration/impl/SharedPointer.hpp"
#include "logging/Logger.hpp"

namespace precice {
namespace acceleration {

class AitkenAcceleration : public Acceleration {
public:
  AitkenAcceleration(
      double                  initialRelaxationFactor,
      std::vector<int>        dataIDs,
      impl::PtrPreconditioner preconditioner);

  virtual ~AitkenAcceleration() {}

  virtual std::vector<int> getDataIDs() const
  {
    return _dataIDs;
  }

  virtual void initialize(
      const DataMap &cpldata);

  virtual void performAcceleration(
      const DataMap &cpldata);

  virtual void iterationsConverged(
      const DataMap &cpldata);

private:
  logging::Logger _log{"acceleration::AitkenAcceleration"};

  double _initialRelaxation;

  std::vector<int> _dataIDs;

  double _aitkenFactor;

  int _iterationCounter = 0;

  Eigen::VectorXd _oldResiduals;
  Eigen::VectorXd _values;
  Eigen::VectorXd _oldValues;

  /// Preconditioner for data vector if multiple data sets are used.
  impl::PtrPreconditioner _preconditioner;
};
} // namespace acceleration
} // namespace precice
