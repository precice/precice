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

  virtual std::vector<int> getPrimaryDataIDs() const override final
  {
    return _primaryDataIDs;
  }

  virtual void initialize(
      const DataMap &cpldata) override final;

  virtual void performAcceleration(
      DataMap &cpldata) override final;

  virtual void iterationsConverged(
      const DataMap &cpldata) override final;

private:
  /// @brief Concatenates the data and old data in cplData into two long vectors
  void concatenateCouplingData(
      const DataMap &cplData, const std::vector<DataID> &dataIDs, Eigen::VectorXd &targetValues, Eigen::VectorXd &targetOldValues) const;

  logging::Logger _log{"acceleration::AitkenAcceleration"};

  const double _initialRelaxation;

  const std::vector<int> _primaryDataIDs;

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
