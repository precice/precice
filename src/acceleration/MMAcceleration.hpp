#pragma once

#include <Eigen/Core>
#include <deque>

#include "logging/Logger.hpp"
#include "acceleration/Acceleration.hpp"
#include "acceleration/SharedPointer.hpp"
#include "acceleration/impl/Preconditioner.hpp"
#include "acceleration/impl/SharedPointer.hpp"

namespace precice
{
namespace acceleration
{

/// Base Class for quasi-Newton acceleration schemes
class MMAcceleration : public Acceleration
{
public:
  MMAcceleration(
      PtrAcceleration       coarseModelOptimization,
      int                     maxIterationsUsed,
      int                     timestepsReused,
      int                     filter,
      double                  singularityLimit,
      bool                    estimateJacobian,
      std::vector<int>        fineDataIDs,
      std::vector<int>        coarseDataIDs,
      impl::PtrPreconditioner       preconditioner);

  virtual ~MMAcceleration()
  {
  }

  /// Returns all MM involved fine model data IDs.
  virtual std::vector<int> getDataIDs() const
  {
    return _fineDataIDs;
  }

  /// Returns all MM involved coarse model data IDs.
  std::vector<int> getCoarseDataIDs() const
  {
    return _coarseDataIDs;
  }

  /// Initializes the acceleration.
  virtual void initialize(DataMap &cplData);

  /**
   * @brief Performs one acceleration step.
   *
   * Has to be called after every implicit coupling iteration.
   */
  virtual void performAcceleration(DataMap &cplData);

  /**
   * @brief Marks a iteration sequence as converged.
   *
   * Since convergence measurements are done outside the acceleration, this
   * method has to be used to signalize convergence to the acceleration.
   */
  virtual void iterationsConverged(DataMap &cplData);

  /**
   * @brief sets the design specification we want to meet for the objective function,
   *     i. e., we want to solve for argmin_x ||R(x) - q||, with R(x) = H(x) - x
   *     Usually we want to solve for a fixed-point of H, thus solving for argmin_x ||R(x)||
   *     with q=0.
   */
  virtual void setDesignSpecification(Eigen::VectorXd &q);

  /**
   * @brief Returns the design specification for the optimization problem.
   *        Information needed to measure the convergence.
   *        In case of manifold mapping it also returns the design specification
   *        for the surrogate model which is updated in every iteration.
   */
  virtual std::map<int, Eigen::VectorXd> getDesignSpecification(DataMap &cplData);

  /**
   * @brief Sets whether the solver has to evaluate the coarse or the fine model representation
   * steers the coupling scheme and the acceleration.
   */
  virtual void setCoarseModelOptimizationActive(bool *coarseOptActive)
  {
    _isCoarseModelOptimizationActive = coarseOptActive;
  }

  /// Exports the current state of the acceleration to a file.
  virtual void exportState(io::TXTWriter &writer);

  /**
   * @brief Imports the last exported state of the acceleration from file.
   *
   * Is empty at the moment!!!
   */
  virtual void importState(io::TXTReader &reader);

  // delete this:
  virtual int getDeletedColumns();

  /// Indicates whether the given acceleration is based on a multi-level approach
  virtual bool isMultilevelBasedApproach()
  {
    return true;
  }

private:
  logging::Logger _log{"acceleration::MMAcceleration"};

  /// Coarse model optimization method
  PtrAcceleration _coarseModelOptimization;

  /// preconditioner, i.e., scaling operator for the LS system
  impl::PtrPreconditioner _preconditioner;

  /// Maximum number of old data iterations kept.
  int _maxIterationsUsed;

  /// Maximum number of old timesteps (with data values) kept.
  int _timestepsReused;

  // @brief Determines sensitivity when two matrix columns are considered equal.
  //
  // When during the QR decomposition of the V matrix a pivot element smaller
  // than the singularity limit is found, the matrix is considered to be singular
  // and the corresponding (older) iteration is removed.
  double _singularityLimit;

  /**
   * @brief sets the design specification we want to meet for the objective function,
   *     i. e., we want to solve for argmin_x ||R(x) - q||, with R(x) = H(x) - x
   *     Usually we want to solve for a fixed-point of H, thus solving for argmin_x ||R(x)||
   *     with q=0.
   */
  Eigen::VectorXd _designSpecification;
  Eigen::VectorXd _coarseModel_designSpecification;

  /**
   * @brief Sets whether the solver has to evaluate the coarse or the fine model representation
   * steers the coupling scheme and the acceleration.
   */
  bool *_isCoarseModelOptimizationActive = nullptr;

  /// @brief Data IDs of data to be involved in the MM algorithm.
  std::vector<int> _fineDataIDs;
  std::vector<int> _coarseDataIDs;
  std::vector<int> _dataIDs;

  /// @brief Data IDs of data not involved in MM coefficient computation.
  std::vector<int> _secondaryDataIDs;

  /** @brief Scales data by fixed value.
   *
   * When more than one data is used to compute the PP update, all data should have similar
   * magnitude, in order to be similarly important in the least-squares solution.
   */
  std::map<int, double> _scalings;

  // @brief Indicates the first iteration, where constant relaxation is used.
  bool _firstIteration = true;

  // @brief Indicates the first time step, where constant relaxation is used
  //        later, we replace the constant relaxation by a qN-update from last time step.
  bool _firstTimeStep = true;

  /** @brief Indicates whether the Jacobian is stored explicitly (multi-vector method) or
   *         if a matrix-vector update is used without explicit representation of the Jacobian.
   */
  bool _estimateJacobian;

  /// current iteration residuals of fine model coupling data
  Eigen::VectorXd _fineResiduals;

  /// current iteration residuals of coarse model coupling data
  Eigen::VectorXd _coarseResiduals;

  /// difference between solver input and output of fine model from last time step
  Eigen::VectorXd _fineOldResiduals;

  /// difference between solver input and output of coarse model from last time step
  Eigen::VectorXd _coarseOldResiduals;

  /// Temporary used in performAcceleration(). output fine model
  Eigen::VectorXd _outputFineModel;

  /// Temporary used in performAcceleration().
  //Eigen::VectorXd _scaledOldValues;

  /// Temporary used in performAcceleration(), output coarse model
  Eigen::VectorXd _outputCoarseModel;

  /// Temporary used in performAcceleration(), input for model evaluation
  Eigen::VectorXd _input_Xstar;

  /// Temporary used in performAcceleration().
  //Eigen::VectorXd _coarseScaledOldValues;

  /// Stores residual deltas for the fine model response
  Eigen::MatrixXd _matrixF;

  /// Stores residual deltas for the coarse model response
  Eigen::MatrixXd _matrixC;

  /// The pseudo inverse of the manifold mapping matrix, only stored and updated if _estimateJacobian is set to true.
  Eigen::MatrixXd _MMMappingMatrix;
  Eigen::MatrixXd _MMMappingMatrix_prev;

  /** @brief Indices (of columns in F, C matrices) of 1st iterations of timesteps.
   *
   * When old timesteps are reused (_timestepsReused > 0), the indices of the
   * first iteration of each timestep needs to be stored, such that, e.g., all
   * iterations of the last timestep, or one specific iteration that leads to
   * a singular matrix in the QR decomposition can be removed and tracked.
   */
  std::deque<int> _matrixCols;

  /** @brief only needed for the parallel master-slave mode. stores the local dimensions,
   *        i.e., the offsets in _invJacobian for all processors
   */
  std::vector<int> _dimOffsets;

  int _iterCoarseModelOpt = 0;
  int _maxIterCoarseModelOpt;

  /// only debugging info, remove this:
  int its = 0, tSteps = 0;
  int deletedColumns = 0;

  /** @brief filter method that is used to maintain good conditioning of the least-squares system
   *		Either of two types: QR1FILTER or QR2Filter
   */
  int _filter;

  bool _notConvergedWithinMaxIter = false;

  /** @brief: computes number of cols in least squares system, i.e, number of cols in
   * 		  _matrixV, _matrixW, _qrV, etc..
   *		  This is necessary only for master-slave mode, when some procs do not have
   *		  any nodes on the coupling interface. In this case, the matrices are not
   * 		  constructed and we have no information about the number of cols. This info
   * 		  is needed for master-slave communication.
   * 		  Number of its =! _cols in general.
   */
  int getLSSystemCols();
  int getLSSystemRows();

  /// updates the V, W matrices (as well as the matrices for the secondary data)
  void updateDifferenceMatrices(
      DataMap &cplData);

  /** @brief registers the new solution x_k+1 (x_star) from the coarse model optimization
   *         problem as new input data for the fine model evaluation step. This has to be done
   *         in each iteration that performs the coarse model optimization.
   */
  void registerSolutionCoarseModelOptimization(DataMap &cplData);

  /** @brief: computes/updates the design specification for the coarse model optimization problem
   *     	   i. e., q_k = c(x_k) - T_k * (f(x_k) - q), q = 0 is the fine model design specification
   */
  void computeCoarseModelDesignSpecifiaction();

  /// computes the quasi-Newton update using the specified pp scheme (MVQN, IQNILS)
  void computeQNUpdate(DataMap &cplData, Eigen::VectorXd &xUpdate);

  /// Removes one iteration from V,W matrices and adapts _matrixCols.
  void removeMatrixColumn(int columnIndex);

  /// concatenates all coupling data involved in the QN system in a single vector
  void concatenateCouplingData(DataMap &cplData);

  /// Indicates whether the design specification has been set and is active or not
  bool isSet(Eigen::VectorXd &designSpec);
};
}
} // namespace precice, acceleration
