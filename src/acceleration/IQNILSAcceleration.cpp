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
  //initialize x_tildes for secondary data
  for (int id : _secondaryDataIDs) {
    precice::time::Storage localCopy = cplData.at(id)->timeStepsStorage();
    _secondaryOldXTildesW.insert(std::pair<int, precice::time::Storage>(id, localCopy));
  }
  // do common QN acceleration initialization
  BaseQNAcceleration::initialize(cplData);
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
      addSecondaryWaveforms(cplData);
    }
  }

  // call the base method for common update of V, W matrices
  BaseQNAcceleration::updateDifferenceMatrices(cplData);
}

void IQNILSAcceleration::computeUnderrelaxationSecondaryData(
    const DataMap &cplData)
{

  //Store x_tildes for secondary data
  _secondaryOldXTildesW.clear();
  for (int id : _secondaryDataIDs) {
    precice::time::Storage localCopy = cplData.at(id)->timeStepsStorage();
    _secondaryOldXTildesW.insert(std::pair<int, precice::time::Storage>(id, localCopy));
  }

  // Perform underrelaxation with initial relaxation factor for secondary data can use the waveform variant for both cases
  for (int id : _secondaryDataIDs) {
    PtrCouplingData data = cplData.at(id);

    for (auto &stample : data->stamples()) {
      auto values    = stample.sample.values;
      auto oldValues = data->getPreviousValuesAtTime(stample.timestamp); // IMPORTANT DETAIL: The interpolation that we use for resampling does not necessarily have to be the same interpolation as the interpolation the user accesses via read-data. (But probably it is easier to just use the same)
      data->values() = values * _initialRelaxation;
      data->values() += oldValues * (1 - _initialRelaxation);
      // Apply relaxation to all timesteps and store it in the current waveform
      data->setSampleAtTime(stample.timestamp, data->sample());
    }
  }
}

void IQNILSAcceleration::computeQNUpdate(const DataMap &cplData)
{

  //Store x_tildes for secondary data
  _secondaryOldXTildesW.clear();
  for (int id : _secondaryDataIDs) {
    precice::time::Storage localCopy = cplData.at(id)->timeStepsStorage();
    _secondaryOldXTildesW.insert(std::pair<int, precice::time::Storage>(id, localCopy));
  }

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

  /**
       * perform QN-Update step for the waveform iteration
       * This is equivalent QN acceleration without waveforms, since only the last sample is updated in that case
       */

  // If the previous time window converged within one single iteration, nothing was added
  // to the LS system matrices and they need to be restored from the backup at time T-2
  if (not _firstTimeWindow && (getLSSystemCols() < 1) && (_timeWindowsReused == 0) && not _forceInitialRelaxation) {
    PRECICE_DEBUG("   Last time window converged after one iteration. Need to restore the secondaryMatricesW from backup.");
    _secondaryWaveformW = _secondaryWaveformWBackup;
  }

  for (int id : _dataIDs) {

    std::vector<precice::time::Storage> Wlist = _waveformW[id];
    //skip the first sample since it always contains the initial data that never changes
    for (auto &stample : cplData.at(id)->stamples()) {
      cplData.at(id)->values() = stample.sample.values;
      double timestamp         = stample.timestamp;
      for (int i = 0; i < c.size(); i++) {
        cplData.at(id)->values() += Wlist[i].sample(timestamp) * c[i];
      }

      //  if the updates resulted in Nan values
      if ((cplData.at(id)->values().array() != cplData.at(id)->values().array()).any()) {
        PRECICE_ERROR("The quasi-Newton update contains NaN values. This means that the quasi-Newton acceleration failed to converge. "
                      "When writing your own adapter this could indicate that you give wrong information to preCICE, such as identical "
                      "data in succeeding iterations. Or you do not properly save and reload checkpoints. "
                      "If you give the correct data this could also mean that the coupled problem is too hard to solve. Try to use a QR "
                      "filter or increase its threshold (larger epsilon).");
      }
      cplData.at(id)->setSampleAtTime(timestamp, cplData.at(id)->sample());
    }
  }

  //Perform QN acceleration for the whole waveform iteration for the secondary ids
  for (int id : _secondaryDataIDs) {

    std::vector<precice::time::Storage> Wlist = _secondaryWaveformW[id];

    //skip the first sample since it always contains the initial data that never changes
    for (auto &stample : cplData.at(id)->stamples()) {

      double timestamp = stample.timestamp;

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
    _secondaryWaveformW.clear();
  }
}

void IQNILSAcceleration::specializedIterationsConverged(
    const DataMap &cplData)
{
  PRECICE_TRACE();
  if (_matrixCols.empty()) {
    PRECICE_WARN("The IQN matrix has no columns.");
  } else {
    if (_matrixCols.front() == 0) { // Did only one iteration
      _matrixCols.pop_front();
    }
    if (_timeWindowsReused == 0) {
      if (_forceInitialRelaxation) {
        _secondaryWaveformW.clear();
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
        PRECICE_ASSERT(_secondaryWaveformW[id].size() > toRemove, _secondaryWaveformW[id].size(), toRemove, id);
        for (int i = 0; i < toRemove; i++) {
          _secondaryWaveformW[id].erase(_secondaryWaveformW[id].end());
        }
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
    _secondaryWaveformW[id].erase(_secondaryWaveformW[id].begin() + columnIndex);
  }
  BaseQNAcceleration::removeMatrixColumn(columnIndex);
}

void IQNILSAcceleration::addSecondaryWaveforms(const DataMap &cplData)
{

  bool columnLimitReached = getLSSystemCols() == _maxIterationsUsed;
  bool overdetermined     = getLSSystemCols() <= getLSSystemRows();

  for (int id : _secondaryDataIDs) {

    // Store x_tildes for secondary data
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

} // namespace precice::acceleration
