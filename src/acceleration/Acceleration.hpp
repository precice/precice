#pragma once

#include <Eigen/Core>
#include <map>
#include <vector>

#include "cplscheme/BaseCouplingScheme.hpp"
#include "cplscheme/SharedPointer.hpp"

namespace precice
{
namespace io
{
class TXTWriter;
class TXTReader;
}
}

namespace precice
{
namespace acceleration
{

class Acceleration
{
public:
  static const int NOFILTER      = 0;
  static const int QR1FILTER     = 1;
  static const int QR1FILTER_ABS = 2;
  static const int QR2FILTER     = 3;
  static const int PODFILTER     = 4;

  /// Map from data ID to data values.
  using DataMap   = std::map<int, cplscheme::PtrCouplingData>;
  using ValuesMap = std::map<int, Eigen::VectorXd>;

  virtual ~Acceleration() {}

  virtual std::vector<int> getDataIDs() const = 0;

  virtual void initialize(DataMap &cpldata) = 0;

  virtual void performAcceleration(DataMap &cpldata) = 0;

  virtual void iterationsConverged(DataMap &cpldata) = 0;

  /**
   * @brief sets the design specification we want to meet for the objective function,
   *        i. e., we want to solve for argmin_x ||R(x) - q||, with R(x) = H(x) - x
   *        Usually we want to solve for a fixed-point of H, thus solving for argmin_x ||R(x)||
   *        with q=0.
   */
  virtual void setDesignSpecification(Eigen::VectorXd &q) = 0;

  /**
   * @brief Returns the design specification for the optimization problem.
   *        Information needed to measure the convergence.
   *        In case of manifold mapping it also returns the design specification
   *        for the surrogate model which is updated in every iteration.
   */
  virtual ValuesMap getDesignSpecification(DataMap &cplData) = 0;

  /**
   * @brief Sets whether the solver has to evaluate the coarse or the fine model representation
   *        steers the coupling scheme and the acceleration. Only needed for multilevel based PPs.
   */
  virtual void setCoarseModelOptimizationActive(bool *coarseOptimizationActive){};

  virtual void exportState(io::TXTWriter &writer) {}

  virtual void importState(io::TXTReader &reader) {}

  /**
   * @brief performs one optimization step of the optimization problem
   *        x_k = argmin_x||f(x_k) - q_k)
   *        with the design specification q_k and the model response f(x_k)
   */
  virtual void optimize(DataMap &cplData, Eigen::VectorXd &q)
  {
    setDesignSpecification(q);
    performAcceleration(cplData);
  };

  virtual int getDeletedColumns()
  {
    return 0;
  }

  /// Indicates whether the given acceleration is based on a multi-level approach
  virtual bool isMultilevelBasedApproach()
  {
    return false;
  }
};
}
} // namespace precice, acceleration
