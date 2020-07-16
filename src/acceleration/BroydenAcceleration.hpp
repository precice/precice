#pragma once

#include <Eigen/Core>
#include <iosfwd>
#include <vector>
#include "acceleration/Acceleration.hpp"
#include "acceleration/BaseQNAcceleration.hpp"
#include "acceleration/impl/SharedPointer.hpp"

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace acceleration {

/**
 * @brief Multi vector quasi-Newton update scheme 
 *
 * Performs a multi vector quasi-Newton to accelerate the convergence of implicit coupling
 * iterations. A multi Broyden update, together with the reuse of the approximate inverse 
 * Jacobian from the old time step are used to approximate the inverse Jacobian. After every
 * coupling iteration, the data values used are enhanced by the new coupling iterates.
 *
 * If more coupling data is present than used to compute the MVQN acceleration,
 * this data is relaxed using the same linear combination as computed for the
 * MVQN-related data. The data is called "secondary" henceforth and additional
 * old value and data matrices are needed for it.
 */
class BroydenAcceleration : public BaseQNAcceleration {
public:
  /**
   * @brief Constructor.
   */
  BroydenAcceleration(
      double                  initialRelaxation,
      bool                    forceInitialRelaxation,
      int                     maxIterationsUsed,
      int                     timestepsReused,
      int                     filter,
      double                  singularityLimit,
      std::vector<int>        dataIDs,
      impl::PtrPreconditioner preconditioner);

  /**
    * @brief Destructor, empty.
    */
  virtual ~BroydenAcceleration() {}

  /**
    * @brief Initializes the acceleration.
    */
  virtual void initialize(DataMap &cplData);

  /**
    * @brief Marks a iteration sequence as converged.
    *
    * called by the iterationsConverged() method in the BaseQNAcceleration class
    * handles the acceleration sepcific action after the convergence of one iteration
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
  int _currentColumns = 0;

  // @brief comptes the MVQN update using QR decomposition of V,
  //        furthermore it updates the inverse of the system jacobian
  virtual void computeQNUpdate(DataMap &cplData, Eigen::VectorXd &xUpdate);

  // @brief updates the V, W matrices (as well as the matrices for the secondary data)
  virtual void updateDifferenceMatrices(DataMap &cplData);

  // @brief computes underrelaxation for the secondary data
  virtual void computeUnderrelaxationSecondaryData(DataMap &cplData);
  //void computeNewtonFactorsQRDecomposition(DataMap& cplData, Eigen::VectorXd& update);
};
} // namespace acceleration
} // namespace precice
