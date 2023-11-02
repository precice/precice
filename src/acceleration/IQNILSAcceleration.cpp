#include "acceleration/IQNILSAcceleration.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <deque>
#include <memory>
#include <utility>

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
#include "utils/IntraComm.hpp"
#include "utils/assertion.hpp"

//#include "utils/NumericalCompare.hpp"

using precice::cplscheme::PtrCouplingData;

namespace precice::acceleration {

IQNILSAcceleration::IQNILSAcceleration(
    double                  initialRelaxation,
    bool                    forceInitialRelaxation,
    int                     maxIterationsUsed,
    int                     pastTimeWindowsReused,
    int                     filter,
    double                  singularityLimit,
    std::vector<int>        dataIDs,
    impl::PtrPreconditioner preconditioner)
    : BaseQNAcceleration(initialRelaxation, forceInitialRelaxation, maxIterationsUsed, pastTimeWindowsReused,
                         filter, singularityLimit, std::move(dataIDs), std::move(preconditioner))
{
}

void IQNILSAcceleration::initialize(
    const DataMap &cplData)
{
  // do common QN acceleration initialization
  BaseQNAcceleration::initialize(cplData);

  // Fetch secondary data IDs, to be relaxed with same coefficients from IQN-ILS
  for (const DataMap::value_type &pair : cplData) {
    if (not utils::contained(pair.first, _dataIDs)) {
      int secondaryEntries = pair.second->getSize();
      utils::append(_secondaryOldXTildes[pair.first], Eigen::VectorXd(Eigen::VectorXd::Zero(secondaryEntries)));
    }
  }
}

void IQNILSAcceleration::updateDifferenceMatrices(
    const DataMap &cplData)
{
  // Compute residuals of secondary data
  for (int id : _secondaryDataIDs) {
    Eigen::VectorXd &secResiduals = _secondaryResiduals[id];
    PtrCouplingData  data         = cplData.at(id);
    PRECICE_ASSERT(secResiduals.size() == data->getSize(),
                   secResiduals.size(), data->getSize());
    secResiduals = data->values();
    secResiduals -= data->previousIteration();
  }

  if (_firstIteration && (_firstTimeWindow || _forceInitialRelaxation)) {
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
        PRECICE_ASSERT(secW.rows() == cplData.at(id)->getSize(), secW.rows(), cplData.at(id)->getSize());
        secW.col(0) = cplData.at(id)->values();
        secW.col(0) -= _secondaryOldXTildes[id];

        std::vector<precice::time::Storage> &vec       = _secondaryWaveformW[id];
        precice::time::Storage               localCopy = cplData.at(id)->timeStepsStorage();

        if (columnLimitReached || overdetermined) {
          vec.erase(vec.end());
        }

        for (auto stample : localCopy.stamples()) {
          stample.sample.values -= _secondaryOldXTildesW.at(id).sample(stample.timestamp);
          localCopy.setSampleAtTime(stample.timestamp, stample.sample);
        }
        vec.insert(vec.begin(), localCopy);
      }
    }

    // Store x_tildes for secondary data
    for (int id : _secondaryDataIDs) {
      PRECICE_ASSERT(_secondaryOldXTildes[id].size() == cplData.at(id)->getSize(),
                     _secondaryOldXTildes[id].size(), cplData.at(id)->getSize());
      _secondaryOldXTildes[id] = cplData.at(id)->values();

      _secondaryOldXTildesW.clear();
      for (int id : _dataIDs) {
        precice::time::Storage localCopy = cplData.at(id)->timeStepsStorage();
        _secondaryOldXTildesW.insert(std::pair<int, precice::time::Storage>(id, localCopy));
      }
    }
  }

  // call the base method for common update of V, W matrices
  BaseQNAcceleration::updateDifferenceMatrices(cplData);
}

void IQNILSAcceleration::computeUnderrelaxationSecondaryData(
    const DataMap &cplData)
{
  //Store x_tildes for secondary data
  for (int id : _secondaryDataIDs) {
    PRECICE_ASSERT(_secondaryOldXTildes.at(id).size() == cplData.at(id)->getSize(),
                   _secondaryOldXTildes.at(id).size(), cplData.at(id)->getSize());
    _secondaryOldXTildes[id] = cplData.at(id)->values();
  }

  // Perform underrelaxation with initial relaxation factor for secondary data
  for (int id : _secondaryDataIDs) {
    PtrCouplingData  data   = cplData.at(id);
    Eigen::VectorXd &values = data->values();
    values *= _initialRelaxation; // new * omg
    Eigen::VectorXd &secResiduals = _secondaryResiduals[id];
    secResiduals                  = data->previousIteration(); // old
    secResiduals *= 1.0 - _initialRelaxation;                  // (1-omg) * old
    values += secResiduals;                                    // (1-omg) * old + new * omg
  }

  for (const DataMap::value_type &pair : cplData) {
    const auto couplingData = pair.second;

    for (auto &stample : couplingData->stamples()) {
      auto values            = stample.sample.values;
      auto oldValues         = couplingData->getPreviousValuesAtTime(stample.timestamp); // IMPORTANT DETAIL: The interpolation that we use for resampling does not necessarily have to be the same interpolation as the interpolation the user accesses via read-data. (But probably it is easier to just use the same)
      couplingData->values() = values * _initialRelaxation;
      couplingData->values() += oldValues * (1 - _initialRelaxation);
      // Apply relaxation to all timesteps and store it in the current waveform
      couplingData->setSampleAtTime(stample.timestamp, couplingData->sample());
    }
  }
}

void IQNILSAcceleration::computeQNUpdate(const DataMap &cplData, Eigen::VectorXd &xUpdate)
{
  PRECICE_TRACE();
  PRECICE_DEBUG("   Compute Newton factors");

  // Calculate QR decomposition of matrix V and solve Rc = -Qr
  Eigen::VectorXd c;

  // for procs with no vertices,
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
  utils::append(c, Eigen::VectorXd(Eigen::VectorXd::Zero(_local_b.size())));

  // compute rhs Q^T*res in parallel
  if (!utils::IntraComm::isParallel()) {
    PRECICE_ASSERT(Q.cols() == getLSSystemCols(), Q.cols(), getLSSystemCols());
    // back substitution
    c = R.triangularView<Eigen::Upper>().solve<Eigen::OnTheLeft>(_local_b);
  } else {
    PRECICE_ASSERT(utils::IntraComm::getCommunication() != nullptr);
    PRECICE_ASSERT(utils::IntraComm::getCommunication()->isConnected());
    if (_hasNodesOnInterface) {
      PRECICE_ASSERT(Q.cols() == getLSSystemCols(), Q.cols(), getLSSystemCols());
    }
    PRECICE_ASSERT(_local_b.size() == getLSSystemCols(), _local_b.size(), getLSSystemCols());

    if (utils::IntraComm::isPrimary()) {
      PRECICE_ASSERT(_global_b.size() == 0, _global_b.size());
    }
    utils::append(_global_b, Eigen::VectorXd(Eigen::VectorXd::Zero(_local_b.size())));

    // do a reduce operation to sum up all the _local_b vectors
    utils::IntraComm::reduceSum(_local_b, _global_b);

    // back substitution R*c = b only on the primary rank
    if (utils::IntraComm::isPrimary()) {
      c = R.triangularView<Eigen::Upper>().solve<Eigen::OnTheLeft>(_global_b);
    }

    // broadcast coefficients c to all secondary ranks
    utils::IntraComm::broadcast(c);
  }

  PRECICE_DEBUG("   Apply Newton factors");
  // compute x updates from W and coefficients c, i.e, xUpdate = c*W
  xUpdate = _matrixW * c;

  /**
     *  perform QN-Update step for the secondary Data
     */

  // If the previous time window converged within one single iteration, nothing was added
  // to the LS system matrices and they need to be restored from the backup at time T-2
  if (not _firstTimeWindow && (getLSSystemCols() < 1) && (_timeWindowsReused == 0) && not _forceInitialRelaxation) {
    PRECICE_DEBUG("   Last time window converged after one iteration. Need to restore the secondaryMatricesW from backup.");
    _secondaryMatricesW = _secondaryMatricesWBackup;
    _secondaryWaveformW = _secondaryWaveformWBackup;
  }

  // Perform QN relaxation for secondary data
  for (int id : _secondaryDataIDs) {
    PtrCouplingData data   = cplData.at(id);
    auto &          values = data->values();
    PRECICE_ASSERT(_secondaryMatricesW[id].cols() == c.size(), _secondaryMatricesW[id].cols(), c.size());
    values = _secondaryMatricesW[id] * c;
    PRECICE_ASSERT(data->getSize() == data->getPreviousIterationSize(), data->getSize(), data->getPreviousIterationSize());
    values += data->previousIteration();
    PRECICE_ASSERT(data->getSize() == _secondaryResiduals[id].size(), data->getSize(), _secondaryResiduals[id].size());
    values += _secondaryResiduals[id];
  }

  // Perform QN acceleration for the whole waveform iteration
  for (int id : _dataIDs) {

    std::vector<precice::time::Storage> Wlist = _waveformW[id];

    for (auto &stample : cplData.at(id)->stamples()) {

      cplData.at(id)->values() = stample.sample.values;
      double timestamp         = stample.timestamp;
      for (int i = 0; i < c.size(); i++) {
        cplData.at(id)->values() += Wlist[i].sample(timestamp) * c[i];
      }
      cplData.at(id)->setSampleAtTime(timestamp, cplData.at(id)->sample());
    }
  }

  // Perform QN acceleration for the whole waveform iteration for the secondary ids
  for (int id : _secondaryDataIDs) {

    std::vector<precice::time::Storage> Wlist = _secondaryWaveformW[id];
    for (auto &stample : cplData.at(id)->stamples()) {

      cplData.at(id)->values() = stample.sample.values;
      double timestamp         = stample.timestamp;

      for (int i = 0; i < c.size(); i++) {
        cplData.at(id)->values() += Wlist[i].sample(timestamp) * c[i];
      }
      cplData.at(id)->setSampleAtTime(timestamp, cplData.at(id)->sample());
    }
  }

  // pending deletion: delete old secondaryMatricesW
  if (_firstIteration && _timeWindowsReused == 0 && not _forceInitialRelaxation) {
    // save current secondaryMatrix data in case the coupling for the next time window will terminate
    // after the first iteration (no new data, i.e., V = W = 0)
    if (getLSSystemCols() > 0) {
      _secondaryWaveformWBackup = _secondaryWaveformW;
    }

    for (int id : _secondaryDataIDs) {
      _secondaryMatricesW[id].resize(0, 0);
      _secondaryWaveformW.clear();
    }
  }
}

void IQNILSAcceleration::specializedIterationsConverged(
    const DataMap &cplData)
{
  PRECICE_TRACE();
  if (_matrixCols.front() == 0) { // Did only one iteration
    _matrixCols.pop_front();
  }

  if (_timeWindowsReused == 0) {
    if (_forceInitialRelaxation) {
      for (int id : _secondaryDataIDs) {
        _secondaryMatricesW[id].resize(0, 0);
      }
    } else {
      /**
       * pending deletion (after first iteration of next time window
       * Using the matrices from the old time window for the first iteration
       * is better than doing underrelaxation as first iteration of every time window
       */
    }
  } else if (static_cast<int>(_matrixCols.size()) > _timeWindowsReused) {
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
    _secondaryWaveformW[id].erase(_secondaryWaveformW[id].begin() + columnIndex);
  }

  BaseQNAcceleration::removeMatrixColumn(columnIndex);
}
} // namespace precice::acceleration
