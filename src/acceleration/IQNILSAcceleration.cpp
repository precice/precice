#include "acceleration/IQNILSAcceleration.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <deque>
#include <memory>
#include "acceleration/impl/Preconditioner.hpp"
#include "acceleration/impl/QRFactorization.hpp"
#include "acceleration/impl/SharedPointer.hpp"
#include "com/Communication.hpp"
#include "com/SharedPointer.hpp"
#include "cplscheme/CouplingData.hpp"
#include "cplscheme/SharedPointer.hpp"
#include "logging/LogMacros.hpp"
#include "utils/EigenHelperFunctions.hpp"
#include "utils/Helpers.hpp"
#include "utils/MasterSlave.hpp"
#include "utils/assertion.hpp"

//#include "utils/NumericalCompare.hpp"

using precice::cplscheme::PtrCouplingData;

namespace precice {
namespace acceleration {

IQNILSAcceleration::IQNILSAcceleration(
    double                  initialRelaxation,
    bool                    forceInitialRelaxation,
    int                     maxIterationsUsed,
    int                     timestepsReused,
    int                     filter,
    double                  singularityLimit,
    std::vector<int>        dataIDs,
    impl::PtrPreconditioner preconditioner)
    : BaseQNAcceleration(initialRelaxation, forceInitialRelaxation, maxIterationsUsed, timestepsReused,
                         filter, singularityLimit, dataIDs, preconditioner)
{
}

void IQNILSAcceleration::initialize(
    DataMap &cplData)
{
  // do common QN acceleration initialization
  BaseQNAcceleration::initialize(cplData);

  // Fetch secondary data IDs, to be relaxed with same coefficients from IQN-ILS
  for (DataMap::value_type &pair : cplData) {
    if (not utils::contained(pair.first, _dataIDs)) {
      int secondaryEntries = pair.second->values().size();
      utils::append(_secondaryOldXTildes[pair.first], (Eigen::VectorXd) Eigen::VectorXd::Zero(secondaryEntries));
    }
  }
}

void IQNILSAcceleration::updateDifferenceMatrices(
    DataMap &cplData)
{
  // Compute residuals of secondary data
  for (int id : _secondaryDataIDs) {
    Eigen::VectorXd &secResiduals = _secondaryResiduals[id];
    PtrCouplingData  data         = cplData[id];
    PRECICE_ASSERT(secResiduals.size() == data->values().size(),
                   secResiduals.size(), data->values().size());
    secResiduals = data->values();
    secResiduals -= data->oldValues.col(0);
  }

  if (_firstIteration && (_firstTimeStep || _forceInitialRelaxation)) {
    // constant relaxation: for secondary data called from base class
  } else {
    if (not _firstIteration) {
      bool columnLimitReached = getLSSystemCols() == _maxIterationsUsed;
      bool overdetermined     = getLSSystemCols() <= getLSSystemRows();
      if (not columnLimitReached && overdetermined) {

        // Append column for secondary W matrices
        for (int id : _secondaryDataIDs) {
          utils::appendFront(_secondaryMatricesW[id], _secondaryResiduals[id]);
        }
      } else {
        // Shift column for secondary W matrices
        for (int id : _secondaryDataIDs) {
          utils::shiftSetFirst(_secondaryMatricesW[id], _secondaryResiduals[id]);
        }
      }

      // Compute delta_x_tilde for secondary data
      for (int id : _secondaryDataIDs) {
        Eigen::MatrixXd &secW = _secondaryMatricesW[id];
        PRECICE_ASSERT(secW.rows() == cplData[id]->values().size(), secW.rows(), cplData[id]->values().size());
        secW.col(0) = cplData[id]->values();
        secW.col(0) -= _secondaryOldXTildes[id];
      }
    }

    // Store x_tildes for secondary data
    for (int id : _secondaryDataIDs) {
      PRECICE_ASSERT(_secondaryOldXTildes[id].size() == cplData[id]->values().size(),
                     _secondaryOldXTildes[id].size(), cplData[id]->values().size());
      _secondaryOldXTildes[id] = cplData[id]->values();
    }
  }

  // call the base method for common update of V, W matrices
  BaseQNAcceleration::updateDifferenceMatrices(cplData);
}

void IQNILSAcceleration::computeUnderrelaxationSecondaryData(
    DataMap &cplData)
{
  //Store x_tildes for secondary data
  for (int id : _secondaryDataIDs) {
    PRECICE_ASSERT(_secondaryOldXTildes[id].size() == cplData[id]->values().size(),
                   _secondaryOldXTildes[id].size(), cplData[id]->values().size());
    _secondaryOldXTildes[id] = cplData[id]->values();
  }

  // Perform underrelaxation with initial relaxation factor for secondary data
  for (int id : _secondaryDataIDs) {
    PtrCouplingData  data   = cplData[id];
    Eigen::VectorXd &values = data->values();
    values *= _initialRelaxation; // new * omg
    Eigen::VectorXd &secResiduals = _secondaryResiduals[id];
    secResiduals                  = data->oldValues.col(0); // old
    secResiduals *= 1.0 - _initialRelaxation;               // (1-omg) * old
    values += secResiduals;                                 // (1-omg) * old + new * omg
  }
}

void IQNILSAcceleration::computeQNUpdate(Acceleration::DataMap &cplData, Eigen::VectorXd &xUpdate)
{
  PRECICE_TRACE();
  PRECICE_DEBUG("   Compute Newton factors");

  // Calculate QR decomposition of matrix V and solve Rc = -Qr
  Eigen::VectorXd c;

  // for master-slave mode and procs with no vertices,
  // qrV.cols() = getLSSystemCols() and _qrV.rows() = 0
  auto Q = _qrV.matrixQ();
  auto R = _qrV.matrixR();

  if (!_hasNodesOnInterface) {
    PRECICE_ASSERT(_qrV.cols() == getLSSystemCols(), _qrV.cols(), getLSSystemCols());
    PRECICE_ASSERT(_qrV.rows() == 0, _qrV.rows());
    PRECICE_ASSERT(Q.size() == 0, Q.size());
  }

  Eigen::VectorXd _local_b = Eigen::VectorXd::Zero(_qrV.cols());
  Eigen::VectorXd _global_b;

  // need to scale the residual to compensate for the scaling in c = R^-1 * Q^T * P^-1 * residual'
  // it is also possible to apply the inverse scaling weights from the right to the vector c
  _preconditioner->apply(_residuals);
  _local_b = Q.transpose() * _residuals;
  _preconditioner->revert(_residuals);
  _local_b *= -1.0; // = -Qr

  PRECICE_ASSERT(c.size() == 0, c.size());
  // reserve memory for c
  utils::append(c, (Eigen::VectorXd) Eigen::VectorXd::Zero(_local_b.size()));

  // compute rhs Q^T*res in parallel
  if (not utils::MasterSlave::isMaster() && not utils::MasterSlave::isSlave()) {
    PRECICE_ASSERT(Q.cols() == getLSSystemCols(), Q.cols(), getLSSystemCols());
    // back substitution
    c = R.triangularView<Eigen::Upper>().solve<Eigen::OnTheLeft>(_local_b);
  } else {
    PRECICE_ASSERT(utils::MasterSlave::_communication.get() != nullptr);
    PRECICE_ASSERT(utils::MasterSlave::_communication->isConnected());
    if (_hasNodesOnInterface) {
      PRECICE_ASSERT(Q.cols() == getLSSystemCols(), Q.cols(), getLSSystemCols());
    }
    PRECICE_ASSERT(_local_b.size() == getLSSystemCols(), _local_b.size(), getLSSystemCols());

    if (utils::MasterSlave::isMaster()) {
      PRECICE_ASSERT(_global_b.size() == 0, _global_b.size());
    }
    utils::append(_global_b, (Eigen::VectorXd) Eigen::VectorXd::Zero(_local_b.size()));

    // do a reduce operation to sum up all the _local_b vectors
    utils::MasterSlave::reduceSum(_local_b.data(), _global_b.data(), _local_b.size()); // size = getLSSystemCols() = _local_b.size()

    // back substitution R*c = b only in master node
    if (utils::MasterSlave::isMaster())
      c = R.triangularView<Eigen::Upper>().solve<Eigen::OnTheLeft>(_global_b);

    // broadcast coefficients c to all slaves
    utils::MasterSlave::broadcast(c.data(), c.size());
  }

  PRECICE_DEBUG("   Apply Newton factors");
  // compute x updates from W and coefficients c, i.e, xUpdate = c*W
  xUpdate = _matrixW * c;

  //PRECICE_DEBUG("c = " << c);

  /**
     *  perform QN-Update step for the secondary Data
     */

  // If the previous time step converged within one single iteration, nothing was added
  // to the LS system matrices and they need to be restored from the backup at time T-2
  if (not _firstTimeStep && (getLSSystemCols() < 1) && (_timestepsReused == 0) && not _forceInitialRelaxation) {
    PRECICE_DEBUG("   Last time step converged after one iteration. Need to restore the secondaryMatricesW from backup.");
    _secondaryMatricesW = _secondaryMatricesWBackup;
  }

  // Perform QN relaxation for secondary data
  for (int id : _secondaryDataIDs) {
    PtrCouplingData data   = cplData[id];
    auto &          values = data->values();
    PRECICE_ASSERT(_secondaryMatricesW[id].cols() == c.size(), _secondaryMatricesW[id].cols(), c.size());
    values = _secondaryMatricesW[id] * c;
    PRECICE_ASSERT(values.size() == data->oldValues.col(0).size(), values.size(), data->oldValues.col(0).size());
    values += data->oldValues.col(0);
    PRECICE_ASSERT(values.size() == _secondaryResiduals[id].size(), values.size(), _secondaryResiduals[id].size());
    values += _secondaryResiduals[id];
  }

  // pending deletion: delete old secondaryMatricesW
  if (_firstIteration && _timestepsReused == 0 && not _forceInitialRelaxation) {
    // save current secondaryMatrix data in case the coupling for the next time step will terminate
    // after the first iteration (no new data, i.e., V = W = 0)
    if (getLSSystemCols() > 0) {
      _secondaryMatricesWBackup = _secondaryMatricesW;
    }
    for (int id : _secondaryDataIDs) {
      _secondaryMatricesW[id].resize(0, 0);
    }
  }
}

void IQNILSAcceleration::specializedIterationsConverged(
    DataMap &cplData)
{
  PRECICE_TRACE();
  if (_matrixCols.front() == 0) { // Did only one iteration
    _matrixCols.pop_front();
  }

  if (_timestepsReused == 0) {
    if (_forceInitialRelaxation) {
      for (int id : _secondaryDataIDs) {
        _secondaryMatricesW[id].resize(0, 0);
      }
    } else {
      /**
       * pending deletion (after first iteration of next time step
       * Using the matrices from the old time step for the first iteration
       * is better than doing underrelaxation as first iteration of every time step
       */
    }
  } else if ((int) _matrixCols.size() > _timestepsReused) {
    int toRemove = _matrixCols.back();
    for (int id : _secondaryDataIDs) {
      Eigen::MatrixXd &secW = _secondaryMatricesW[id];
      PRECICE_ASSERT(secW.cols() > toRemove, secW, toRemove, id);
      for (int i = 0; i < toRemove; i++) {
        utils::removeColumnFromMatrix(secW, secW.cols() - 1);
      }
    }
  }
}

void IQNILSAcceleration::removeMatrixColumn(
    int columnIndex)
{
  PRECICE_ASSERT(_matrixV.cols() > 1);
  // remove column from secondary Data Matrix W
  for (int id : _secondaryDataIDs) {
    utils::removeColumnFromMatrix(_secondaryMatricesW[id], columnIndex);
  }

  BaseQNAcceleration::removeMatrixColumn(columnIndex);
}
} // namespace acceleration
} // namespace precice
