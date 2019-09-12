#pragma once

#include <map>

#include "logging/Logger.hpp"
#include "acceleration/Acceleration.hpp"

namespace precice
{
namespace acceleration
{

class ConstantRelaxationAcceleration : public Acceleration
{
public:
  ConstantRelaxationAcceleration(
      double           relaxation,
      std::vector<int> dataIDs);

  virtual ~ConstantRelaxationAcceleration() {}

  virtual std::vector<int> getDataIDs() const
  {
    return _dataIDs;
  }

  virtual void setDesignSpecification(
      Eigen::VectorXd &q);

  virtual std::map<int, Eigen::VectorXd> getDesignSpecification(DataMap &cplData);

  virtual void initialize(DataMap &cplData);

  virtual void performAcceleration(DataMap &cplData);

  virtual void iterationsConverged(DataMap &cplData)
  {
  }

private:
  logging::Logger _log{"acceleration::ConstantRelaxationAcceleration"};

  double _relaxation;

  std::vector<int> _dataIDs;

  Eigen::VectorXd _designSpecification;
};
}
} // namespace precice, acceleration
