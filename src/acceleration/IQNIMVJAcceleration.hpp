#ifndef PRECICE_NO_MPI

#pragma once

#include <Eigen/Core>
#include <deque>
#include <vector>
#include "acceleration/Acceleration.hpp"
#include "acceleration/BaseQNAcceleration.hpp"
#include "acceleration/impl/ParallelMatrixOperations.hpp"
#include "acceleration/impl/SVDFactorization.hpp"
#include "acceleration/impl/SharedPointer.hpp"
#include "com/SharedPointer.hpp"

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice::acceleration {

/**
 * @brief Multi vector quasi-Newton update scheme
 *
 * Performs a multi vector quasi-Newton to accelerate the convergence of implicit coupling
 * iterations. A multi Broyden update, together with the reuse of the approximate inverse
 * Jacobian from the old time window are used to approximate the inverse Jacobian. After every
 * coupling iteration, the data values used are enhanced by the new coupling iterates.
 *
 * If more coupling data is present than used to compute the IQNIMVJ acceleration,
 * this data is relaxed using the same linear combination as computed for the
 * IQNIMVJ-related data. The data is called "secondary" henceforth and additional
 * old value and data matrices are needed for it.
 */
class IQNIMVJAcceleration : public BaseQNAcceleration {
public:
  static const int NO_RESTART = 0;
  static const int RS_ZERO    = 1;
  static const int RS_LS      = 2;
  static const int RS_SVD     = 3;
  static const int RS_SLIDE   = 4;

  /**
   * @brief Constructor.
   */
  IQNIMVJAcceleration(
      double                         initialRelaxation,
      bool                           forceInitialRelaxation,
      int                            maxIterationsUsed,
      int                            pastTimeWindowsReused,
      int                            filter,
      double                         singularityLimit,
      std::vector<int>               dataIDs,
      const impl::PtrPreconditioner &preconditioner,
      bool                           alwaysBuildJacobian,
      int                            imvjRestartType,
      int                            chunkSize,
      int                            RSLSreusedTimeWindows,
      double                         RSSVDtruncationEps,
      bool                           reducedTimeGrid);

  /**
   * @brief Destructor, empty.
   */
  ~IQNIMVJAcceleration() override;

  /**
   * @brief Marks a iteration sequence as converged.
   *
   * called by the iterationsConverged() method in the BaseQNAcceleration class
   * handles the acceleration specific action after the convergence of one iteration
   */
  void specializedIterationsConverged(const DataMap &cplData) override;

private:
  /// @brief stores the approximation of the inverse Jacobian of the system at current time window.
  Eigen::MatrixXd _invJacobian;

  /// @brief stores the approximation of the inverse Jacobian from the previous time window.
  Eigen::MatrixXd _oldInvJacobian;

  /// @brief stores the sub result (W-J_prev*V) for the current iteration
  Eigen::MatrixXd _Wtil;

  /// @brief stores all Wtil matrices within the current chunk of the imvj restart mode, disabled if _imvjRestart = false.
  std::vector<Eigen::MatrixXd> _WtilChunk;

  /// @brief stores all pseudo inverses within the current chunk of the imvj restart mode, disabled if _imvjRestart = false.
  std::vector<Eigen::MatrixXd> _pseudoInverseChunk;

  /// @brief stores columns from previous  #_RSLSreusedTimeWindows time windows if RS-LS restart-mode is active
  Eigen::MatrixXd _matrixV_RSLS;

  /// @brief stores columns from previous  #_RSLSreusedTimeWindows time windows if RS-LS restart-mode is active
  Eigen::MatrixXd _matrixW_RSLS;

  /// @brief number of cols per time window
  std::deque<int> _matrixCols_RSLS;

  /// @brief encapsulates matrix-matrix and matrix-vector multiplications for serial and parallel execution
  impl::PtrParMatrixOps _parMatrixOps;

  /// @brief holds and maintains a truncated SVD decomposition of the Jacobian matrix
  impl::SVDFactorization _svdJ;

  /** @brief If true, the less efficient method to compute the quasi-Newton update is used,
   *   that explicitly builds the Jacobian in each iteration. If set to false this is only done
   *   in the very last iteration and the update is computed based on MATVEC products.
   */
  bool _alwaysBuildJacobian;

  /** @brief: Indicates the type of the imvj restart-mode:
   *  - NO_RESTART: imvj is run on normal mode which builds the Jacobian explicitly
   *  - RS-ZERO:    imvj is run in restart-mode. After M time windows all stored matrices are dropped
   *  - RS-LS:      imvj in restart-mode. After M time windows restart with LS approximation for initial Jacobian
   *  - RS-SVD:     imvj in restart mode. After M time windows, update of an truncated SVD of the Jacobian.
   */
  const int _imvjRestartType;

  /** @brief: If true, the imvj method is used with the restart chunk based approach that avoids
   *  to explicitly build and store the Jacobian. If false, the Jacobian is stored and build, however,
   *  no truncation of information is present.
   */
  const bool _imvjRestart;

  /// @brief: Number of time windows between restarts for the imvj method in restart mode
  int _chunkSize;

  /// @brief: Number of reused time windows at restart if restart-mode = RS-LS
  int _RSLSreusedTimeWindows;

  /// @brief tracks the number of restarts of IMVJ
  int _nbRestarts;

  // DEBUG
  // std::fstream _info2;
  double _avgRank;

  /** @brief: computes the IQNIMVJ update using QR decomposition of V,
   *        furthermore it updates the inverse of the system jacobian
   */
  void computeQNUpdate(Eigen::VectorXd &xUpdate) override;

  /// @brief: updates the V, W matrices (as well as the matrices for the secondary data)
  void updateDifferenceMatrices(const DataMap &cplData) override;

  /** @brief: computes the quasi-Newton update vector based on the matrices V and W using a QR
   *  decomposition of V. The decomposition is not re-computed en-block in every iteration
   *  but updated so that the new added column in V is incorporated in the decomposition.
   *
   *  This method rebuilds the Jacobian matrix and the matrix W_til in each iteration
   *  which is not necessary and thus inefficient.
   */
  void computeNewtonUpdate(Eigen::VectorXd &update);

  /** @brief: computes the quasi-Newton update vector based on the same numerics as above.
   *  However, it exploits the fact that the matrix W_til can be updated according to V and W
   *  via the formula W_til.col(j) = W.col(j) - J_inv * V.col(j).
   *  Then, pure matrix-vector products are sufficient to compute the update within one iteration, i.e.,
   *  (1) x1 := J_prev*(-res) (2) y := Z(-res) (3) xUp := W_til*y + x1
   *  The Jacobian matrix only needs to be set up in the very last iteration of one time window, i.e.
   *  in iterationsConverged.
   */
  void computeNewtonUpdateEfficient(Eigen::VectorXd &update);

  /** @brief: computes the pseudo inverse of V multiplied with V^T, i.e., Z = (V^TV)^-1V^T via QR-dec
   */
  void pseudoInverse(Eigen::MatrixXd &pseudoInverse);

  /** @brief: computes a explicit representation of the Jacobian, i.e., n x n matrix
   */
  void buildJacobian();

  /** @brief: re-computes the matrix _Wtil = ( W - J_prev * V) instead of updating it according to V
   */
  void buildWtil();

  /** @brief: restarts the imvj method, i.e., drops all stored matrices Wtil and Z and computes a
   *  initial guess of the Jacobian based on the given restart strategy:
   *  RS-LS:   Perform a IQN-LS least squares initial guess with _RSLSreusedTimeWindows
   *  RS-SVD:  Update a truncated SVD decomposition of the SVD with rank-1 modifications from Wtil*Z
   *  RS-Zero: Start with zero information, initial guess J = 0.
   */
  void restartIMVJ();

  /// @brief: Removes one iteration from V,W matrices and adapts _matrixCols.
  void removeMatrixColumn(int columnIndex) override;

  /// @copydoc precice::Acceleration::BaseQNAcceleration::specializedInitializeVectorsAndPreconditioner()
  void specializedInitializeVectorsAndPreconditioner(const DataMap &cplData) final;

  /// @brief: Removes one column form the V_RSLS and W_RSLS matrices and adapts _matrixCols_RSLS
  void removeMatrixColumnRSLS(int columnINdex);
};
} // namespace precice::acceleration

#endif /* PRECICE_NO_MPI */
