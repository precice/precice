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

namespace precice::acceleration {

class AitkenAcceleration : public Acceleration {
public:
  AitkenAcceleration(
      double                  initialRelaxationFactor,
      std::vector<int>        dataIDs,
      impl::PtrPreconditioner preconditioner);

  ~AitkenAcceleration() override = default;

  std::vector<int> getPrimaryDataIDs() const final
  {
    return _primaryDataIDs;
  }

  void initialize(
      const DataMap &cpldata) final;

  void performAcceleration(
      DataMap &cpldata,
      double   windowStart,
      double   windowEnd) final;

  void iterationsConverged(
      const DataMap &cpldata, double windowStart) final;

private:
  /// @brief Concatenates the data and old data in cplData into two long vectors
  void concatenateCouplingData(
      const DataMap &cplData, const std::vector<DataID> &dataIDs, Eigen::VectorXd &targetValues, Eigen::VectorXd &targetOldValues, double windowStart, double windowEnd) const;

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
} // namespace precice::acceleration
