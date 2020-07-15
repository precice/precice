#pragma once

#include <Eigen/Core>
#include <algorithm>
#include <deque>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include "acceleration/Acceleration.hpp"
#include "acceleration/impl/QRFactorization.hpp"
#include "acceleration/impl/SharedPointer.hpp"
#include "logging/Logger.hpp"

/* ****************************************************************************
 *
 * A few comments concerning the present design choice.
 *
 * All the functions from the base class BAseQNAcceleration are specialized in
 * the sub classes as needed. This is done vi overwriting the base functions in
 * the specialized sub classes and calling the respective base function after
 * performing the specialized stuff (in order to perform the common, generalized
 * computations i.e. handling of V,W matrices etc.)
 * However, for the performAcceleration Method we decided (for better readability)
 * to have this method only in the base class, while introducing a function
 * performPPSecondaryData that handles all the specialized stuff concerning acceleration
 * processing for the secondary data in the sub classes.
 *
 * Another possibility would have been to introduce a bunch of functions like
 * initializeSpecialized(), removeMatrixColumnSpecialized(),
 * iterationsConvergedSpecialized(), etc
 * and call those function from the base class top down to the sub classes.
 *
 * The third possibility was to separate the approximation of the Jacobian from
 * the common stuff like handling V,W matrices in the acceleration.
 * Here, we have a class QNAcceleration that handles the V,W stuff an d the basic
 * scheme of the QN update. Furthermore we have a base class (or rather interface)
 * JacobianApproximation with sub classes MVQNAPX and IQNAPX that handle all the
 * specialized stuff like Jacobian approximation, handling of secondary data etc.
 * However, this approach is not feasible, as we have to call the function
 * removeMatrixColumn() down in the specialized sub classes MVQNApx and IQNApx.
 * This is not possible as the function works on the V, W matrices that are
 * completely treated by QNAcceleration.
 *
 * ****************************************************************************
 */

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace io {
class TXTReader;
class TXTWriter;
} // namespace io

namespace acceleration {

/**
 * @brief Base Class for quasi-Newton acceleration schemes
 *
 */
class BaseQNAcceleration : public Acceleration {
public:
  BaseQNAcceleration(
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
  virtual ~BaseQNAcceleration()
  {
    // not necessary for user, only for developer, if needed, this should be configurable
    //     if (utils::MasterSlave::isMaster() || (not utils::MasterSlave::isMaster() && not utils::MasterSlave::isSlave())){
    //       _infostream.open("precice-accelerationInfo.log", std::ios_base::out);
    //       _infostream << std::setprecision(16);
    //       _infostream << _infostringstream.str();
    //     }
  }

  /**
    * @brief Returns all IQN involved data IDs.
    */
  virtual std::vector<int> getDataIDs() const
  {
    return _dataIDs;
  }

  /**
    * @brief Initializes the acceleration.
    */
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
    * @brief Exports the current state of the acceleration to a file.
    */
  virtual void exportState(io::TXTWriter &writer);

  /**
    * @brief Imports the last exported state of the acceleration from file.
    *
    * Is empty at the moment!!!
    */
  virtual void importState(io::TXTReader &reader);

  /// how many QN columns were deleted in this timestep
  virtual int getDeletedColumns() const;

  /// how many QN columns were dropped (went out of scope) in this timestep
  virtual int getDroppedColumns() const;

  /** @brief: computes number of cols in least squares system, i.e, number of cols in
    *  _matrixV, _matrixW, _qrV, etc..
    *	 This is necessary only for master-slave mode, when some procs do not have
    *	 any nodes on the coupling interface. In this case, the matrices are not
    *  constructed and we have no information about the number of cols. This info
    *  is needed for master-slave communication. Number of its =! _cols in general.
    */
  virtual int getLSSystemCols() const;

protected:
  logging::Logger _log{"acceleration::BaseQNAcceleration"};

  /// Preconditioner for least-squares system if vectorial system is used.
  impl::PtrPreconditioner _preconditioner;

  /// Constant relaxation factor used for first iteration.
  double _initialRelaxation;

  /// Maximum number of old data iterations kept.
  int _maxIterationsUsed;

  /// Maximum number of old timesteps (with data values) kept.
  int _timestepsReused;

  /// Data IDs of data to be involved in the IQN algorithm.
  std::vector<int> _dataIDs;

  /// Data IDs of data not involved in IQN coefficient computation.
  std::vector<int> _secondaryDataIDs;

  /// Indicates the first iteration, where constant relaxation is used.
  bool _firstIteration = true;

  /* @brief Indicates the first time step, where constant relaxation is used
    *        later, we replace the constant relaxation by a qN-update from last time step.
    */
  bool _firstTimeStep = true;

  /*
    * @brief If in master-slave mode: True if this process has nodes at the coupling interface
    */
  bool _hasNodesOnInterface = true;

  /* @brief If true, the QN-scheme always performs a underrelaxation in the first iteration of
    *        a new time step. Otherwise, the LS system from the previous time step is used in the
    *        first iteration.
    */
  bool _forceInitialRelaxation;

  /** @brief If true, the LS system has been modified (reset or recomputed) in such a way, that mere
    *         updating of matrices _Wtil, Q, R etc.. is not feasible any more and need to be recomputed.
    */
  bool _resetLS = false;

  /// @brief Solver output from last iteration.
  Eigen::VectorXd _oldXTilde;

  /// @brief Current iteration residuals of IQN data. Temporary.
  Eigen::VectorXd _residuals;

  /// @brief Current iteration residuals of secondary data.
  std::map<int, Eigen::VectorXd> _secondaryResiduals;

  /// @brief Stores residual deltas.
  Eigen::MatrixXd _matrixV;

  /// @brief Stores x tilde deltas, where x tilde are values computed by solvers.
  Eigen::MatrixXd _matrixW;

  /// @brief Stores the current QR decomposition ov _matrixV, can be updated via deletion/insertion of columns
  impl::QRFactorization _qrV;

  /** @brief filter method that is used to maintain good conditioning of the least-squares system
    *        Either of two types: QR1FILTER or QR2Filter
    */
  int _filter;

  /** @brief Determines sensitivity when two matrix columns are considered equal.
    *
    * When during the QR decomposition of the V matrix a pivot element smaller
    * than the singularity limit is found, the matrix is considered to be singular
    * and the corresponding (older) iteration is removed.
    */
  double _singularityLimit;

  /** @brief Indices (of columns in W, V matrices) of 1st iterations of timesteps.
    *
    * When old timesteps are reused (_timestepsReused > 0), the indices of the
    * first iteration of each timestep needs to be stored, such that, e.g., all
    * iterations of the last timestep, or one specific iteration that leads to
    * a singular matrix in the QR decomposition can be removed and tracked.
    */
  std::deque<int> _matrixCols;

  /** @brief only needed for the parallel master-slave mode. stores the local dimensions,
    *  i.e., the offsets in _invJacobian for all processors
    */
  std::vector<int> _dimOffsets;

  /// @brief write some debug/acceleration info to file
  std::ostringstream _infostringstream;
  std::fstream       _infostream;

  int getLSSystemRows();

  /**
     * @brief Marks a iteration sequence as converged.
     *
     * called by the iterationsConverged() method in the BaseQNAcceleration class
     * handles the acceleration sepcific action after the convergence of one iteration
     */
  virtual void specializedIterationsConverged(DataMap &cplData) = 0;

  /// Updates the V, W matrices (as well as the matrices for the secondary data)
  virtual void updateDifferenceMatrices(DataMap &cplData);

  /// Concatenates all coupling data involved in the QN system in a single vector
  virtual void concatenateCouplingData(DataMap &cplData);

  /// Splits up QN system vector back into the coupling data
  virtual void splitCouplingData(DataMap &cplData);

  /// Applies the filter method for the least-squares system, defined in the configuration
  virtual void applyFilter();

  /// Computes underrelaxation for the secondary data
  virtual void computeUnderrelaxationSecondaryData(DataMap &cplData) = 0;

  /// Computes the quasi-Newton update using the specified pp scheme (MVQN, IQNILS)
  virtual void computeQNUpdate(DataMap &cplData, Eigen::VectorXd &xUpdate) = 0;

  /// Removes one iteration from V,W matrices and adapts _matrixCols.
  virtual void removeMatrixColumn(int columnIndex);

  /// Wwrites info to the _infostream (also in parallel)
  void writeInfo(std::string s, bool allProcs = false);

  int its = 0, tSteps = 0;

private:
  /// @brief Concatenation of all coupling data involved in the QN system.
  Eigen::VectorXd _values;

  /// @brief Concatenation of all (old) coupling data involved in the QN system.
  Eigen::VectorXd _oldValues;

  /// @brief Difference between solver input and output from last timestep
  Eigen::VectorXd _oldResiduals;

  /** @brief backup of the V,W and matrixCols data structures. Needed for the skipping of
   *  initial relaxation, if previous time step converged within one iteration i.e., V and W
   *  are empty -- in this case restore V and W with time step t-2.
   */
  Eigen::MatrixXd _matrixVBackup;
  Eigen::MatrixXd _matrixWBackup;
  std::deque<int> _matrixColsBackup;

  /// Number of filtered out columns in this time window
  int _nbDelCols = 0;

  /// Number of dropped columns in this time window (old time window out of scope)
  int _nbDropCols = 0;
};
} // namespace acceleration
} // namespace precice
