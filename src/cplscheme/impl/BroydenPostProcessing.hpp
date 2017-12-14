#pragma once

#include "BaseQNPostProcessing.hpp"

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice
{
namespace cplscheme
{
namespace impl
{

/**
 * @brief Multi vector quasi-Newton update scheme 
 *
 * Performs a multi vector quasi-Newton to accelerate the convergence of implicit coupling
 * iterations. A multi Broyden update, together with the reuse of the approximate inverse 
 * Jacobian from the old time step are used to approximate the inverse Jacobian. After every
 * coupling iteration, the data values used are enhanced by the new coupling iterates.
 *
 * If more coupling data is present than used to compute the MVQN post-processing,
 * this data is relaxed using the same linear combination as computed for the
 * MVQN-related data. The data is called "secondary" henceforth and additional
 * old value and data matrices are needed for it.
 */
class BroydenPostProcessing : public BaseQNPostProcessing
{
public:
  /**
   * @brief Constructor.
   */
  BroydenPostProcessing(
      double            initialRelaxation,
      bool              forceInitialRelaxation,
      int               maxIterationsUsed,
      int               timestepsReused,
      int               filter,
      double            singularityLimit,
      std::vector<int>  dataIDs,
      PtrPreconditioner preconditioner);

  /**
    * @brief Destructor, empty.
    */
  virtual ~BroydenPostProcessing() {}

  /**
    * @brief Initializes the post-processing.
    */
  virtual void initialize(DataMap &cplData);

  /**
    * @brief Marks a iteration sequence as converged.
    *
    * called by the iterationsConverged() method in the BaseQNPostProcessing class
    * handles the postprocessing sepcific action after the convergence of one iteration
    */
  virtual void specializedIterationsConverged(DataMap &cplData);

private:
  // remove this ofter debugging, not useful
  // ---------------------------------------
  std::fstream f;
  //----------------------------------------

  // @brief stores the approximation of the inverse Jacobian of the system at current time step.
  Eigen::MatrixXd _invJacobian;
  Eigen::MatrixXd _oldInvJacobian;

  int _maxColumns;
  int _currentColumns;

  // @brief comptes the MVQN update using QR decomposition of V,
  //        furthermore it updates the inverse of the system jacobian
  virtual void computeQNUpdate(DataMap &cplData, Eigen::VectorXd &xUpdate);

  // @brief updates the V, W matrices (as well as the matrices for the secondary data)
  virtual void updateDifferenceMatrices(DataMap &cplData);

  // @brief computes underrelaxation for the secondary data
  virtual void computeUnderrelaxationSecondaryData(DataMap &cplData);
  //void computeNewtonFactorsQRDecomposition(DataMap& cplData, Eigen::VectorXd& update);
};
}
}
} // namespace precice, cplscheme, impl
