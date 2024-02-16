#ifndef PRECICE_NO_MPI

#include <Eigen/Core>
#include <algorithm>
#include <fstream>
#include <map>
#include <memory>
#include <string>
#include <utility>

#include <boost/range/adaptor/map.hpp>
#include "acceleration/IQNIMVJAcceleration.hpp"
#include "acceleration/impl/ParallelMatrixOperations.hpp"
#include "acceleration/impl/Preconditioner.hpp"
#include "acceleration/impl/QRFactorization.hpp"
#include "com/Communication.hpp"
#include "com/MPIPortsCommunication.hpp"
#include "cplscheme/CouplingData.hpp"
#include "cplscheme/SharedPointer.hpp"
#include "logging/LogMacros.hpp"
#include "precice/impl/Types.hpp"
#include "profiling/Event.hpp"
#include "utils/EigenHelperFunctions.hpp"
#include "utils/IntraComm.hpp"
#include "utils/assertion.hpp"

using precice::cplscheme::PtrCouplingData;

namespace precice::acceleration {

// ==================================================================================
IQNIMVJAcceleration::IQNIMVJAcceleration(
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
    double                         RSSVDtruncationEps)
    : BaseQNAcceleration(initialRelaxation, forceInitialRelaxation, maxIterationsUsed, pastTimeWindowsReused,
                         filter, singularityLimit, std::move(dataIDs), preconditioner),
      //  _secondaryOldXTildes(),
      _invJacobian(),
      _oldInvJacobian(),
      _Wtil(),
      _WtilChunk(),
      _pseudoInverseChunk(),
      _matrixV_RSLS(),
      _matrixW_RSLS(),
      _matrixCols_RSLS(),
      _parMatrixOps(nullptr),
      _svdJ(RSSVDtruncationEps, preconditioner),
      _alwaysBuildJacobian(alwaysBuildJacobian),
      _imvjRestartType(imvjRestartType),
      _imvjRestart(false),
      _chunkSize(chunkSize),
      _RSLSreusedTimeWindows(RSLSreusedTimeWindows),
      _nbRestarts(0),
      _avgRank(0)
{
}

// ==================================================================================
IQNIMVJAcceleration::~IQNIMVJAcceleration() = default;

// ==================================================================================
void IQNIMVJAcceleration::initialize(
    const DataMap &cplData)
{
  PRECICE_TRACE();

  // do common QN acceleration initialization
  BaseQNAcceleration::initialize(cplData);

  if (_imvjRestartType > 0)
    _imvjRestart = true;

  // initialize parallel matrix-matrix operation module
  _parMatrixOps = std::make_shared<impl::ParallelMatrixOperations>();
  _parMatrixOps->initialize(not _imvjRestart);
  _svdJ.initialize(_parMatrixOps, getLSSystemRows());

  int entries  = _residuals.size();
  int global_n = 0;

  if (!utils::IntraComm::isParallel()) {
    global_n = entries;
  } else {
    global_n = _dimOffsets.back();
  }

  if (not _imvjRestart) {
    // only need memory for Jacobain if not in restart mode
    _invJacobian    = Eigen::MatrixXd::Zero(global_n, entries);
    _oldInvJacobian = Eigen::MatrixXd::Zero(global_n, entries);
  }
  // initialize V, W matrices for the LS restart
  if (_imvjRestartType == RS_LS) {
    _matrixCols_RSLS.push_front(0);
    _matrixV_RSLS = Eigen::MatrixXd::Zero(entries, 0);
    _matrixW_RSLS = Eigen::MatrixXd::Zero(entries, 0);
  }
  _Wtil = Eigen::MatrixXd::Zero(entries, 0);

  if (utils::IntraComm::isPrimary() || !utils::IntraComm::isParallel()) {
    _infostringstream << " IMVJ restart mode: " << _imvjRestart << "\n chunk size: " << _chunkSize << "\n trunc eps: " << _svdJ.getThreshold() << "\n R_RS: " << _RSLSreusedTimeWindows << "\n--------\n"
                      << '\n';
  }
  for (const auto &data : cplData | boost::adaptors::map_values) {
    if (data->exchangeSubsteps()) {
      PRECICE_ERROR("Quasi-Newton IMVJ acceleration does not yet support using data from all substeps. Please set substeps=\"false\" in the exchange tag of data \"{}\" or switch to Quasi-Newton ILS acceleration.", data->getDataName());
    }
  }
}
// ==================================================================================
void IQNIMVJAcceleration::computeUnderrelaxationSecondaryData(
    const DataMap &cplData)
{
  // Perform underrelaxation with initial relaxation factor for secondary data
  for (int id : _secondaryDataIDs) {
    PtrCouplingData  data   = cplData.at(id);
    Eigen::VectorXd &values = data->values();
    values *= _initialRelaxation; // new * omg
    Eigen::VectorXd &secResiduals = _secondaryResiduals[id];
    secResiduals                  = data->previousIteration();
    secResiduals *= 1.0 - _initialRelaxation; // (1-omg) * old
    values += secResiduals;                   // (1-omg) * old + new * omg
  }
}

// ==================================================================================
void IQNIMVJAcceleration::updateDifferenceMatrices(
    const DataMap &cplData)
{
  /**
   *  Matrices and vectors used in this method as well as the result Wtil are
   *  NOT SCALED by the preconditioner.
   */

  PRECICE_TRACE();

  // call the base method for common update of V, W matrices
  // important that base method is called before updating _Wtil
  BaseQNAcceleration::updateDifferenceMatrices(cplData);

  // update _Wtil if the efficient computation of the quasi-Newton update is used
  // or update current Wtil if the restart mode of imvj is used
  if (not _alwaysBuildJacobian || _imvjRestart) {
    if (_firstIteration && (_firstTimeWindow || _forceInitialRelaxation)) {
      // do nothing: constant relaxation
    } else {
      if (not _firstIteration) {
        // Update matrix _Wtil = (W - J_prev*V) with newest information

        Eigen::VectorXd v = _matrixV.col(0);

        // Need to reconstruct the newest w vector from _waveformW
        Eigen::VectorXd w = Eigen::VectorXd::Zero(_matrixV.rows());

        // keep track of how far we have filled the vector
        int pos = 0;
        for (int id : _dataIDs) {
          double          endTime                         = _waveformW[id][0].maxStoredTime();
          Eigen::VectorXd lastSample                      = _waveformW[id][0].sample(endTime);
          w(Eigen::seq(pos, pos + lastSample.size() - 1)) = lastSample;
          pos += lastSample.size();
        }

        // here, we check for _Wtil.cols() as the matrices V, W need to be updated before hand
        // and thus getLSSystemCols() does not yield the correct result.
        bool columnLimitReached = _Wtil.cols() == _maxIterationsUsed;
        bool overdetermined     = _Wtil.cols() <= getLSSystemRows();

        Eigen::VectorXd wtil = Eigen::VectorXd::Zero(_matrixV.rows());

        // add column: Wtil(:,0) = W(:,0) - sum_q [ Wtil^q * ( Z^q * V(:,0)) ]
        //                                         |--- J_prev ---|
        // iterate over all stored Wtil and Z matrices in current chunk
        if (_imvjRestart) {
          for (int i = 0; i < static_cast<int>(_WtilChunk.size()); i++) {
            int colsLSSystemBackThen = _pseudoInverseChunk[i].rows();
            PRECICE_ASSERT(colsLSSystemBackThen == _WtilChunk[i].cols(), colsLSSystemBackThen, _WtilChunk[i].cols());
            Eigen::VectorXd Zv = Eigen::VectorXd::Zero(colsLSSystemBackThen);
            // multiply: Zv := Z^q * V(:,0) of size (m x 1)
            _parMatrixOps->multiply(_pseudoInverseChunk[i], v, Zv, colsLSSystemBackThen, getLSSystemRows(), 1);
            // multiply: Wtil^q * Zv  dimensions: (n x m) * (m x 1), fully local
            wtil += _WtilChunk[i] * Zv;
          }

          // store columns if restart mode = RS-LS
          if (_imvjRestartType == RS_LS) {
            utils::appendFront(_matrixV_RSLS, v);
            utils::appendFront(_matrixW_RSLS, w);
            _matrixCols_RSLS.front()++;
          }

          // imvj without restart is used, but efficient update, i.e. no Jacobian assembly in each iteration
          // add column: Wtil(:,0) = W(:,0) - J_prev * V(:,0)
        } else {
          // compute J_prev * V(0) := wtil the new column in _Wtil of dimension: (n x n) * (n x 1) = (n x 1),
          //                                        parallel: (n_global x n_local) * (n_local x 1) = (n_local x 1)
          _parMatrixOps->multiply(_oldInvJacobian, v, wtil, _dimOffsets, getLSSystemRows(), getLSSystemRows(), 1, false);
        }
        wtil *= -1;
        wtil += w;

        if (not columnLimitReached && overdetermined) {
          utils::appendFront(_Wtil, wtil);
        } else {
          utils::shiftSetFirst(_Wtil, wtil);
        }
      }
    }
  }
}

// ==================================================================================
void IQNIMVJAcceleration::computeQNUpdate(
    const DataMap &cplData)
{

  /**
   * The inverse Jacobian
   *
   *        J_inv = J_inv_n + (W - J_inv_n*V)*(V^T*V)^-1*V^T
   *
   * is computed and the resulting quasi-Newton update is returned.
   * Used matrices (V, W, Wtil, invJacobian, oldINvJacobian) are
   * scaled with the used preconditioner.
   */

  /**
   *  The MVJ- quasi-Newton update
   *  Either compute efficient, omitting to build the Jacobian in each iteration or inefficient.
   *  INVARIANT: All objects, J_inv, J_old_inv, W, V, Wtil, xUpdate, res, etc. are unscaled.
   *             Only the QR-decomposition of V is scaled and thus needs to be unscaled before
   *             using it in multiplications with the other matrices.
   */
  Eigen::VectorXd xUpdate = Eigen::VectorXd::Zero(_residuals.size());

  if (_alwaysBuildJacobian) {
    computeNewtonUpdate(cplData, xUpdate);
  } else {
    computeNewtonUpdateEfficient(cplData, xUpdate);
  }
  if (std::isnan(utils::IntraComm::l2norm(xUpdate))) {
    PRECICE_ERROR("The quasi-Newton update contains NaN values. This means that the quasi-Newton acceleration failed to converge. "
                  "When writing your own adapter this could indicate that you give wrong information to preCICE, such as identical "
                  "data in succeeding iterations. Or you do not properly save and reload checkpoints. "
                  "If you give the correct data this could also mean that the coupled problem is too hard to solve. Try to use a QR "
                  "filter or increase its threshold (larger epsilon).");
  }
  _values += xUpdate;

  //Split the coupling data and update the values in coupling data
  BaseQNAcceleration::splitCouplingData(cplData);
}

// ==================================================================================
void IQNIMVJAcceleration::pseudoInverse(
    Eigen::MatrixXd &pseudoInverse)
{
  PRECICE_TRACE();
  /**
   *   computation of pseudo inverse matrix Z = (V^TV)^-1 * V^T as solution
   *   to the equation R*z = Q^T(i) for all columns i,  via back substitution.
   */
  auto Q = _qrV.matrixQ();
  auto R = _qrV.matrixR();

  PRECICE_ASSERT(pseudoInverse.rows() == _qrV.cols(), pseudoInverse.rows(), _qrV.cols());
  PRECICE_ASSERT(pseudoInverse.cols() == _qrV.rows(), pseudoInverse.cols(), _qrV.rows());

  Eigen::VectorXd yVec(pseudoInverse.rows());

  // assertions for the case of processors with no vertices
  if (!_hasNodesOnInterface) {
    PRECICE_ASSERT(_qrV.cols() == getLSSystemCols(), _qrV.cols(), getLSSystemCols());
    PRECICE_ASSERT(_qrV.rows() == 0, _qrV.rows());
    PRECICE_ASSERT(Q.size() == 0, Q.size());
  }

  // backsubstitution
  for (int i = 0; i < Q.rows(); i++) {
    Eigen::VectorXd Qrow = Q.row(i);
    yVec                 = R.triangularView<Eigen::Upper>().solve<Eigen::OnTheLeft>(Qrow);
    pseudoInverse.col(i) = yVec;
  } // ----------------

  // scale pseudo inverse back Z := Z' * P,
  // Z' is scaled pseudo inverse i.e, Z' = R^-1 * Q^T * P^-1
  _preconditioner->apply(pseudoInverse, true);
  //  e.stop(true);
}

// ==================================================================================
void IQNIMVJAcceleration::buildWtil()
{
  /**
   * PRECONDITION: Assumes that V, W, J_prev are already preconditioned,
   */
  PRECICE_TRACE();

  PRECICE_ASSERT(_matrixV.rows() == _qrV.rows(), _matrixV.rows(), _qrV.rows());
  PRECICE_ASSERT(getLSSystemCols() == _qrV.cols(), getLSSystemCols(), _qrV.cols());

  _Wtil = Eigen::MatrixXd::Zero(_qrV.rows(), _qrV.cols());

  // imvj restart mode: re-compute Wtil: Wtil = W - sum_q [ Wtil^q * (Z^q*V) ]
  //                                                      |--- J_prev ---|
  // iterate over all stored Wtil and Z matrices in current chunk
  if (_imvjRestart) {
    for (int i = 0; i < static_cast<int>(_WtilChunk.size()); i++) {
      int colsLSSystemBackThen = _pseudoInverseChunk[i].rows();
      PRECICE_ASSERT(colsLSSystemBackThen == _WtilChunk[i].cols(), colsLSSystemBackThen, _WtilChunk[i].cols());
      Eigen::MatrixXd ZV = Eigen::MatrixXd::Zero(colsLSSystemBackThen, _qrV.cols());
      // multiply: ZV := Z^q * V of size (m x m) with m=#cols, stored on each proc.
      _parMatrixOps->multiply(_pseudoInverseChunk[i], _matrixV, ZV, colsLSSystemBackThen, getLSSystemRows(), _qrV.cols());
      // multiply: Wtil^q * ZV  dimensions: (n x m) * (m x m), fully local and embarrassingly parallel
      _Wtil += _WtilChunk[i] * ZV;
    }

    // imvj without restart is used, i.e., recompute Wtil: Wtil = W - J_prev * V
  } else {
    // multiply J_prev * V = W_til of dimension: (n x n) * (n x m) = (n x m),
    //                                    parallel:  (n_global x n_local) * (n_local x m) = (n_local x m)
    _parMatrixOps->multiply(_oldInvJacobian, _matrixV, _Wtil, _dimOffsets, getLSSystemRows(), getLSSystemRows(), getLSSystemCols(), false);
  }

  // W_til = (W-J_inv_n*V) = (W-V_tilde)
  _Wtil *= -1.;

  // The W matrix is no longer stored as a matrix in the base class but instead as a map from the DataIDs to the k waveforms in W so we need to assemble W matrix

  for (int i = 0; i < _Wtil.cols(); i++) {
    Eigen::VectorXd colVecW = Eigen::VectorXd::Zero(_matrixV.rows());
    // keep track of how far we have filled the vector
    int pos = 0;
    for (int id : _dataIDs) {

      double          endTime                               = _waveformW[id][i].maxStoredTime();
      Eigen::VectorXd lastSample                            = _waveformW[id][i].sample(endTime);
      colVecW(Eigen::seq(pos, pos + lastSample.size() - 1)) = lastSample;
      pos += lastSample.size();
    }
    _Wtil.col(i) = colVecW - _Wtil.col(i);
  }

  _resetLS = false;
  //  e.stop(true);
}

// ==================================================================================
void IQNIMVJAcceleration::buildJacobian()
{
  PRECICE_TRACE();
  /**      --- compute inverse Jacobian ---
   *
   * J_inv = J_inv_n + (W - J_inv_n*V)*(V^T*V)^-1*V^T
   *
   * assumes that J_prev, V, W are already preconditioned
   */

  /**
   *  (1) computation of pseudo inverse Z = (V^TV)^-1 * V^T
   */
  Eigen::MatrixXd Z(_qrV.cols(), _qrV.rows());
  pseudoInverse(Z);

  /**
   *  (2) Multiply J_prev * V =: W_tilde
   */
  PRECICE_ASSERT(_matrixV.rows() == _qrV.rows(), _matrixV.rows(), _qrV.rows());
  PRECICE_ASSERT(getLSSystemCols() == _qrV.cols(), getLSSystemCols(), _qrV.cols());
  if (_resetLS) {
    buildWtil();
    PRECICE_WARN(" ATTENTION, in buildJacobian call for buildWtill() - this should not be the case except the coupling did only one iteration");
  }

  /** (3) compute invJacobian = W_til*Z
   *
   *  where Z = (V^T*V)^-1*V^T via QR-dec and back-substitution       dimension: (n x n) * (n x m) = (n x m),
   *  and W_til = (W - J_inv_n*V)                                     parallel:  (n_global x n_local) * (n_local x m) = (n_local x m)
   */
  _parMatrixOps->multiply(_Wtil, Z, _invJacobian, _dimOffsets, getLSSystemRows(), getLSSystemCols(), getLSSystemRows());
  // --------

  // update Jacobian
  _invJacobian = _invJacobian + _oldInvJacobian;
}

// ==================================================================================
void IQNIMVJAcceleration::computeNewtonUpdateEfficient(
    const DataMap &  cplData,
    Eigen::VectorXd &xUpdate)
{
  PRECICE_TRACE();

  /**      --- update inverse Jacobian efficient, ---
   *   If normal mode is used:
   *   Do not recompute W_til in every iteration and do not build
   *   the entire Jacobian matrix. This is only necessary if the coupling
   *   iteration has converged, namely in the last iteration.
   *
   *   If restart-mode is used:
   *   The Jacobian is never build. Store matrices Wtil^q and Z^q for the last M time windows.
   *   After M time windows, a restart algorithm is performed basedon the restart-mode type, either
   *   Least-Squares restart (IQN-ILS like) or maintaining of a updated truncated SVD decomposition
   *   of the SVD.
   *
   * J_inv = J_inv_n + (W - J_inv_n*V)*(V^T*V)^-1*V^T
   *
   * ASSUMPTION: All objects are scaled with the active preconditioner
   */

  /**
   *  (1) computation of pseudo inverse Z = (V^TV)^-1 * V^T
   */
  Eigen::MatrixXd Z(_qrV.cols(), _qrV.rows());
  pseudoInverse(Z);

  /**
   *  (2) Construction of _Wtil = (W - J_prev * V), should be already present due to updated computation
   */
  PRECICE_ASSERT(_matrixV.rows() == _qrV.rows(), _matrixV.rows(), _qrV.rows());
  PRECICE_ASSERT(getLSSystemCols() == _qrV.cols(), getLSSystemCols(), _qrV.cols());

  // rebuild matrix Wtil if V changes substantially.
  if (_resetLS) {
    buildWtil();
  }

  /**
   *  Avoid computation of Z*Wtil = Jtil \in (n x n). Rather do matrix-vector computations
   *  [ J_prev*(-res) ] + [Wtil * [Z * (-res)] ]
   *  '----- 1 -------'           '----- 2 ----'
   *                      '-------- 3 ---------'
   *
   *  (3) compute r_til = Z*(-residual)   where Z = (V^T*V)^-1*V^T via QR-dec and back-substitution
   *
   *  dimension: (m x n) * (n x 1) = (m x 1),
   *  parallel:  (m x n_local) * (n x 1) = (m x 1)
   */
  Eigen::VectorXd negativeResiduals = -_residuals;
  Eigen::VectorXd r_til             = Eigen::VectorXd::Zero(getLSSystemCols());
  _parMatrixOps->multiply(Z, negativeResiduals, r_til, getLSSystemCols(), getLSSystemRows(), 1); // --------

  /**
   * (4) compute _Wtil * r_til
   *
   * dimension: (n x m) * (m x 1) = (n x 1),
   * parallel:  (n_local x m) * (m x 1) = (n_local x 1)
   *
   * Note: r_til is not distributed but locally stored on each proc (dimension m x 1)
   */
  Eigen::VectorXd xUptmp(_residuals.size());
  xUpdate = Eigen::VectorXd::Zero(_residuals.size());
  xUptmp  = _Wtil * r_til; // local product, result is naturally distributed.

  /**
   *  (5) xUp = J_prev * (-res) + Wtil*Z*(-res)
   *
   *  restart-mode: sum_q { Wtil^q * [ Z^q * (-res) ] },
   *  where r_til = Z^q * (-res) is computed first and then xUp := Wtil^q * r_til
   */
  if (_imvjRestart) {
    for (int i = 0; i < static_cast<int>(_WtilChunk.size()); i++) {
      int colsLSSystemBackThen = _pseudoInverseChunk[i].rows();
      PRECICE_ASSERT(colsLSSystemBackThen == _WtilChunk[i].cols(), colsLSSystemBackThen, _WtilChunk[i].cols());
      r_til = Eigen::VectorXd::Zero(colsLSSystemBackThen);
      // multiply: r_til := Z^q * (-res) of size (m x 1) with m=#cols of LS at that time, result stored on each proc.
      _parMatrixOps->multiply(_pseudoInverseChunk[i], negativeResiduals, r_til, colsLSSystemBackThen, getLSSystemRows(), 1);
      // multiply: Wtil^q * r_til  dimensions: (n x m) * (m x 1), fully local and embarrassingly parallel
      xUpdate += _WtilChunk[i] * r_til;
    }

    // imvj without restart is used, i.e., compute directly J_prev * (-res)
  } else {
    _parMatrixOps->multiply(_oldInvJacobian, negativeResiduals, xUpdate, _dimOffsets, getLSSystemRows(), getLSSystemRows(), 1, false);
    PRECICE_DEBUG("Mult J*V DONE");
  }

  xUpdate += xUptmp;

  // pending deletion: delete Wtil
  if (_firstIteration && _timeWindowsReused == 0 && not _forceInitialRelaxation) {
    _Wtil.conservativeResize(0, 0);
    _resetLS = true;
  }
}

// ==================================================================================
void IQNIMVJAcceleration::computeNewtonUpdate(const DataMap &cplData, Eigen::VectorXd &xUpdate)
{
  PRECICE_TRACE();

  /**      --- update inverse Jacobian ---
   *
   * J_inv = J_inv_n + (W - J_inv_n*V)*(V^T*V)^-1*V^T
   */

  /**  (1) computation of pseudo inverse Z = (V^TV)^-1 * V^T
   */
  Eigen::MatrixXd Z(_qrV.cols(), _qrV.rows());
  pseudoInverse(Z);

  /**  (2) Multiply J_prev * V =: W_tilde
   */
  buildWtil();

  /**  (3) compute invJacobian = W_til*Z
   *
   *  where Z = (V^T*V)^-1*V^T via QR-dec and back-substitution             dimension: (n x n) * (n x m) = (n x m),
   *  and W_til = (W - J_inv_n*V)                                           parallel:  (n_global x n_local) * (n_local x m) = (n_local x m)
   */
  _parMatrixOps->multiply(_Wtil, Z, _invJacobian, _dimOffsets, getLSSystemRows(), getLSSystemCols(), getLSSystemRows()); // --------

  // update Jacobian
  _invJacobian = _invJacobian + _oldInvJacobian;

  /**  (4) solve delta_x = - J_inv * res
   */
  Eigen::VectorXd negativeResiduals = -_residuals;

  // multiply J_inv * (-res) = x_Update of dimension: (n x n) * (n x 1) = (n x 1),
  //                                        parallel: (n_global x n_local) * (n_local x 1) = (n_local x 1)
  _parMatrixOps->multiply(_invJacobian, negativeResiduals, xUpdate, _dimOffsets, getLSSystemRows(), getLSSystemRows(), 1, false); // --------
}

// ==================================================================================
void IQNIMVJAcceleration::restartIMVJ()
{
  PRECICE_TRACE();

  // int used_storage = 0;
  // int theoreticalJ_storage = 2*getLSSystemRows()*_residuals.size() + 3*_residuals.size()*getLSSystemCols() + _residuals.size()*_residuals.size();
  //                ------------ RESTART SVD ------------
  if (_imvjRestartType == IQNIMVJAcceleration::RS_SVD) {

    // we need to compute the updated SVD of the scaled Jacobian matrix
    // |= APPLY PRECONDITIONING  J_prev = Wtil^q, Z^q  ===|
    for (int i = 0; i < static_cast<int>(_WtilChunk.size()); i++) {
      _preconditioner->apply(_WtilChunk[i]);
      _preconditioner->revert(_pseudoInverseChunk[i], true);
    }
    // |===================                            ===|

    Rank rankBefore = _svdJ.isSVDinitialized() ? _svdJ.rank() : 0;

    // if it is the first time window, there is no initial SVD, so take all Wtil, Z matrices
    // otherwise, the first element of each container holds the decomposition of the current
    // truncated SVD, i.e., Wtil^0 = \phi, Z^0 = S\psi^T, this should not be added to the SVD.
    int q = _svdJ.isSVDinitialized() ? 1 : 0;

    // perform M-1 rank-1 updates of the truncated SVD-dec of the Jacobian
    for (; q < static_cast<int>(_WtilChunk.size()); q++) {
      // update SVD, i.e., PSI * SIGMA * PHI^T <-- PSI * SIGMA * PHI^T + Wtil^q * Z^q
      _svdJ.update(_WtilChunk[q], _pseudoInverseChunk[q].transpose());
      //  used_storage += 2*_WtilChunk.size();
    }
    // int m = _WtilChunk[q].cols(), n = _WtilChunk[q].rows();
    // used_storage += 2*rankBefore*m + 4*m*n + 2*m*m + (rankBefore+m)*(rankBefore+m) + 2*n*(rankBefore+m);

    // drop all stored Wtil^q, Z^q matrices
    _WtilChunk.clear();
    _pseudoInverseChunk.clear();

    auto &psi   = _svdJ.matrixPsi();
    auto &sigma = _svdJ.singularValues();
    auto &phi   = _svdJ.matrixPhi();

    // multiply sigma * phi^T, phi is distributed block-row wise, phi^T is distributed block-column wise
    // sigma is stored local on each proc, thus, the multiplication is fully local, no communication.
    // Z = sigma * phi^T
    Eigen::MatrixXd Z(phi.cols(), phi.rows());
    for (int i = 0; i < static_cast<int>(Z.rows()); i++)
      for (int j = 0; j < static_cast<int>(Z.cols()); j++)
        Z(i, j) = phi(j, i) * sigma[i];

    Rank rankAfter = _svdJ.rank();
    int  waste     = _svdJ.getWaste();
    _avgRank += rankAfter;

    // store factorized truncated SVD of J
    _WtilChunk.push_back(psi);
    _pseudoInverseChunk.push_back(Z);

    // |= REVERT PRECONDITIONING  J_prev = Wtil^0, Z^0  ==|
    _preconditioner->revert(_WtilChunk.front());
    _preconditioner->apply(_pseudoInverseChunk.front(), true);
    // |===================                             ==|

    PRECICE_DEBUG("MVJ-RESTART, mode=SVD. Rank of truncated SVD of Jacobian {}, new modes: {}, truncated modes: {} avg rank: {}", rankAfter, rankAfter - rankBefore, waste, _avgRank / _nbRestarts);

    // double percentage = 100.0*used_storage/(double)theoreticalJ_storage;
    if (utils::IntraComm::isPrimary() || !utils::IntraComm::isParallel()) {
      _infostringstream << " - MVJ-RESTART " << _nbRestarts << ", mode= SVD -\n  new modes: " << rankAfter - rankBefore << "\n  rank svd: " << rankAfter << "\n  avg rank: " << _avgRank / _nbRestarts << "\n  truncated modes: " << waste << "\n"
                        << '\n';
    }

    //        ------------ RESTART LEAST SQUARES ------------
  } else if (_imvjRestartType == IQNIMVJAcceleration::RS_LS) {
    // drop all stored Wtil^q, Z^q matrices
    _WtilChunk.clear();
    _pseudoInverseChunk.clear();

    if (_matrixV_RSLS.cols() > 0) {
      // avoid that the system is getting too squared
      while (_matrixV_RSLS.cols() * 2 >= getLSSystemRows()) {
        removeMatrixColumnRSLS(_matrixV_RSLS.cols() - 1);
      }

      // preconditioning
      // V needs to be sclaed to compute the pseudo inverse
      // W only needs to be scaled, as the design requires to store scaled
      // matrices Wtil^0 and Z^0 as initial guess after the restart
      _preconditioner->apply(_matrixV_RSLS);
      _preconditioner->apply(_matrixW_RSLS);

      impl::QRFactorization qr(_filter);
      qr.setGlobalRows(getLSSystemRows());
      // for QR2-filter, the QR-dec is computed in qr-applyFilter()
      if (_filter != Acceleration::QR2FILTER) {
        for (int i = 0; i < static_cast<int>(_matrixV_RSLS.cols()); i++) {
          Eigen::VectorXd v = _matrixV_RSLS.col(i);
          qr.pushBack(v); // same order as matrix V_RSLS
        }
      }

      // apply filter
      if (_filter != Acceleration::NOFILTER) {
        std::vector<int> delIndices(0);
        qr.applyFilter(_singularityLimit, delIndices, _matrixV_RSLS);
        // start with largest index (as V,W matrices are shrunk and shifted
        for (int i = delIndices.size() - 1; i >= 0; i--) {
          removeMatrixColumnRSLS(delIndices[i]);
        }
        PRECICE_ASSERT(_matrixV_RSLS.cols() == qr.cols(), _matrixV_RSLS.cols(), qr.cols());
      }

      /**
       *   computation of pseudo inverse matrix Z = (V^TV)^-1 * V^T as solution
       *   to the equation R*z = Q^T(i) for all columns i,  via back substitution.
       */
      auto            Q = qr.matrixQ();
      auto            R = qr.matrixR();
      Eigen::MatrixXd pseudoInverse(qr.cols(), qr.rows());
      Eigen::VectorXd yVec(pseudoInverse.rows());

      // backsubstitution
      for (int i = 0; i < Q.rows(); i++) {
        Eigen::VectorXd Qrow = Q.row(i);
        yVec                 = R.triangularView<Eigen::Upper>().solve<Eigen::OnTheLeft>(Qrow);
        pseudoInverse.col(i) = yVec;
      }

      // scale pseudo inverse back Z := Z' * P,
      // Z' is scaled pseudo inverse i.e, Z' = R^-1 * Q^T * P^-1
      //_preconditioner->apply(pseudoInverse, true, false);

      // store factorization of least-squares initial guess for Jacobian
      _WtilChunk.push_back(_matrixW_RSLS);
      _pseudoInverseChunk.push_back(pseudoInverse);

      // |= REVERT PRECONDITIONING  J_prev = Wtil^0, Z^0  ==|
      _preconditioner->revert(_WtilChunk.front());
      _preconditioner->apply(_pseudoInverseChunk.front(), true);
      _preconditioner->revert(_matrixW_RSLS);
      _preconditioner->revert(_matrixV_RSLS);
      // |===================                             ==|
    }

    PRECICE_DEBUG("MVJ-RESTART, mode=LS. Restart with {} columns from {} time windows.", _matrixV_RSLS.cols(), _RSLSreusedTimeWindows);
    if (utils::IntraComm::isPrimary() || !utils::IntraComm::isParallel()) {
      _infostringstream << " - MVJ-RESTART" << _nbRestarts << ", mode= LS -\n  used cols: " << _matrixV_RSLS.cols() << "\n  R_RS: " << _RSLSreusedTimeWindows << "\n"
                        << '\n';
    }

    //            ------------ RESTART ZERO ------------
  } else if (_imvjRestartType == IQNIMVJAcceleration::RS_ZERO) {
    // drop all stored Wtil^q, Z^q matrices
    _WtilChunk.clear();
    _pseudoInverseChunk.clear();

    PRECICE_DEBUG("MVJ-RESTART, mode=Zero");

  } else if (_imvjRestartType == IQNIMVJAcceleration::RS_SLIDE) {

    // re-compute Wtil -- compensate for dropping of Wtil_0 and Z_0:
    //                    Wtil_q <-- Wtil_q +  Wtil^0 * (Z^0*V_q)
    for (int i = static_cast<int>(_WtilChunk.size()) - 1; i >= 1; i--) {

      int colsLSSystemBackThen = _pseudoInverseChunk.front().rows();
      PRECICE_ASSERT(colsLSSystemBackThen == _WtilChunk.front().cols(), colsLSSystemBackThen, _WtilChunk.front().cols());
      Eigen::MatrixXd ZV = Eigen::MatrixXd::Zero(colsLSSystemBackThen, _qrV.cols());
      // multiply: ZV := Z^q * V of size (m x m) with m=#cols, stored on each proc.
      _parMatrixOps->multiply(_pseudoInverseChunk.front(), _matrixV, ZV, colsLSSystemBackThen, getLSSystemRows(), _qrV.cols());
      // multiply: Wtil^0 * (Z_0*V)  dimensions: (n x m) * (m x m), fully local and embarrassingly parallel
      Eigen::MatrixXd tmp = Eigen::MatrixXd::Zero(_qrV.rows(), _qrV.cols());
      tmp                 = _WtilChunk.front() * ZV;
      _WtilChunk[i] += tmp;

      // drop oldest pair Wtil_0 and Z_0
      PRECICE_ASSERT(not _WtilChunk.empty());
      PRECICE_ASSERT(not _pseudoInverseChunk.empty())
      _WtilChunk.erase(_WtilChunk.begin());
      _pseudoInverseChunk.erase(_pseudoInverseChunk.begin());
    }

  } else if (_imvjRestartType == IQNIMVJAcceleration::NO_RESTART) {
    PRECICE_ASSERT(false); // should not happen, in this case _imvjRestart=false
  } else {
    PRECICE_ASSERT(false);
  }
}

// ==================================================================================
void IQNIMVJAcceleration::specializedIterationsConverged(
    const DataMap &cplData)
{
  PRECICE_TRACE();

  // truncate V_RSLS and W_RSLS matrices according to _RSLSreusedTimeWindows
  if (_imvjRestartType == RS_LS) {
    if (_matrixCols_RSLS.front() == 0) { // Did only one iteration
      _matrixCols_RSLS.pop_front();
    }
    if (_RSLSreusedTimeWindows == 0) {
      _matrixV_RSLS.resize(0, 0);
      _matrixW_RSLS.resize(0, 0);
      _matrixCols_RSLS.clear();
    } else if (static_cast<int>(_matrixCols_RSLS.size()) > _RSLSreusedTimeWindows) {
      int toRemove = _matrixCols_RSLS.back();
      PRECICE_ASSERT(toRemove > 0, toRemove);
      if (_matrixV_RSLS.size() > 0) {
        PRECICE_ASSERT(_matrixV_RSLS.cols() > toRemove, _matrixV_RSLS.cols(), toRemove);
      }

      // remove columns
      for (int i = 0; i < toRemove; i++) {
        utils::removeColumnFromMatrix(_matrixV_RSLS, _matrixV_RSLS.cols() - 1);
        utils::removeColumnFromMatrix(_matrixW_RSLS, _matrixW_RSLS.cols() - 1);
      }
      _matrixCols_RSLS.pop_back();
    }
    _matrixCols_RSLS.push_front(0);
  }

  //_info2<<'\n';

  // if efficient update of imvj is enabled
  if (not _alwaysBuildJacobian || _imvjRestart) {
    // need to apply the preconditioner, as all data structures are reverted after
    // call to computeQNUpdate. Need to call this before the preconditioner is updated.

    // |= REBUILD QR-dec if needed     ============|
    // apply scaling to V, V' := P * V (only needed to reset the QR-dec of V)
    _preconditioner->apply(_matrixV);

    if (_preconditioner->requireNewQR()) {
      if (not(_filter == Acceleration::QR2FILTER)) { // for QR2 filter, there is no need to do this twice
        _qrV.reset(_matrixV, getLSSystemRows());
      }
      _preconditioner->newQRfulfilled();
    }
    // apply the configured filter to the LS system
    // as it changed in BaseQNAcceleration::iterationsConverged()
    BaseQNAcceleration::applyFilter();
    _preconditioner->revert(_matrixV);
    // |===================          ============|

    //              ------- RESTART/ JACOBIAN ASSEMBLY -------
    if (_imvjRestart) {

      // add the matrices Wtil and Z of the converged configuration to the storage containers
      Eigen::MatrixXd Z(_qrV.cols(), _qrV.rows());
      // compute pseudo inverse using QR factorization and back-substitution
      // also compensates for the scaling of V, i.e.,
      // reverts Z' = R^-1 * Q^T * P^-1 as Z := Z' * P
      pseudoInverse(Z);

      // push back unscaled pseudo Inverse, Wtil is also unscaled.
      // all objects in Wtil chunk and Z chunk are NOT PRECONDITIONED
      _WtilChunk.push_back(_Wtil);
      _pseudoInverseChunk.push_back(Z);

      /**
       *  Restart the IMVJ according to restart type
       */
      if (static_cast<int>(_WtilChunk.size()) >= _chunkSize + 1) {

        // < RESTART >
        _nbRestarts++;
        profiling::Event restartUpdate("IMVJRestart");
        restartIMVJ();
        restartUpdate.stop();
      }

      // only in imvj normal mode with efficient update:
    } else {

      // compute explicit representation of Jacobian
      buildJacobian();
    }

    /** in case of enforced initial relaxation, the matrices are cleared
     *  in case of pastTimeWindowsReused > 0, the columns in _Wtil are outdated, as the Jacobian changes, hence clear
     *  in case of pastTimeWindowsReused == 0 and no initial relaxation, pending deletion in performAcceleration
     */
    if (_timeWindowsReused > 0 || (_timeWindowsReused == 0 && _forceInitialRelaxation)) {
      //_Wtil.conservativeResize(0, 0);
      _resetLS = true;
    }
  }

  // only store Jacobian if imvj is in normal mode, i.e., the Jacobian is build
  if (not _imvjRestart) {
    // store inverse Jacobian from converged time window. NOT SCALED with preconditioner
    _oldInvJacobian = _invJacobian;
  }
}

// ==================================================================================
void IQNIMVJAcceleration::removeMatrixColumn(
    int columnIndex)
{
  PRECICE_TRACE(columnIndex, _matrixV.cols());
  PRECICE_ASSERT(_matrixV.cols() > 1, _matrixV.cols());
  PRECICE_ASSERT(_Wtil.cols() > 1);

  // remove column from matrix _Wtil
  if (not _resetLS && not _alwaysBuildJacobian)
    utils::removeColumnFromMatrix(_Wtil, columnIndex);

  BaseQNAcceleration::removeMatrixColumn(columnIndex);
}

// ==================================================================================
void IQNIMVJAcceleration::removeMatrixColumnRSLS(
    int columnIndex)
{
  PRECICE_TRACE(columnIndex, _matrixV_RSLS.cols());
  PRECICE_ASSERT(_matrixV_RSLS.cols() > 1);

  utils::removeColumnFromMatrix(_matrixV_RSLS, columnIndex);
  utils::removeColumnFromMatrix(_matrixW_RSLS, columnIndex);

  // Reduce column count
  std::deque<int>::iterator iter = _matrixCols_RSLS.begin();
  int                       cols = 0;
  while (iter != _matrixCols_RSLS.end()) {
    cols += *iter;
    if (cols > columnIndex) {
      PRECICE_ASSERT(*iter > 0);
      *iter -= 1;
      if (*iter == 0) {
        _matrixCols_RSLS.erase(iter);
      }
      break;
    }
    iter++;
  }
}

} // namespace precice::acceleration

#endif
