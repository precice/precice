#pragma once

#include <map>
#include "PostProcessing.hpp"
#include "logging/Logger.hpp"

namespace precice
{
namespace cplscheme
{
namespace impl
{

class ConstantRelaxationPostProcessing : public PostProcessing
{
public:
  ConstantRelaxationPostProcessing(
      double           relaxation,
      std::vector<int> dataIDs);

  virtual ~ConstantRelaxationPostProcessing() {}

  virtual std::vector<int> getDataIDs() const
  {
    return _dataIDs;
  }

  virtual void setDesignSpecification(
      Eigen::VectorXd &q);

  virtual std::map<int, Eigen::VectorXd> getDesignSpecification(DataMap &cplData);

  virtual void initialize(DataMap &cplData);

  virtual void performPostProcessing(DataMap &cplData);

  virtual void iterationsConverged(DataMap & /* cplData */)
  {
  }

private:
  logging::Logger _log{"cplscheme::ConstantRelaxationPostProcessing"};

  double _relaxation;

  std::vector<int> _dataIDs;

  Eigen::VectorXd _designSpecification;
};
}
}
} // namespace precice, cplscheme, impl
