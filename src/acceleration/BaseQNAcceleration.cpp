#include "acceleration/BaseQNAcceleration.hpp"
#include <Eigen/Core>
#include <boost/range/adaptor/map.hpp>
#include <cmath>
#include <iomanip>
#include <memory>
#include <utility>

#include "acceleration/impl/Preconditioner.hpp"
#include "acceleration/impl/QRFactorization.hpp"
#include "com/Communication.hpp"
#include "com/SharedPointer.hpp"
#include "cplscheme/CouplingData.hpp"
#include "logging/LogMacros.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "profiling/Event.hpp"
#include "time/TimeGrids.hpp"
#include "utils/EigenHelperFunctions.hpp"
#include "utils/Helpers.hpp"
#include "utils/IntraComm.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace io {
class TXTReader;
class TXTWriter;
} // namespace io
namespace acceleration {

/* ----------------------------------------------------------------------------
 *     Constructor
 * ----------------------------------------------------------------------------
 */
BaseQNAcceleration::BaseQNAcceleration(
    double                  initialRelaxation,
    bool                    forceInitialRelaxation,
    int                     maxIterationsUsed,
    int                     timeWindowsReused,
    int                     filter,
    double                  singularityLimit,
    std::vector<int>        dataIDs,
    std::string             onBoundViolation,
    impl::PtrPreconditioner preconditioner,
    bool                    reducedTimeGrid)
    : _preconditioner(std::move(preconditioner)),
      _initialRelaxation(initialRelaxation),
      _maxIterationsUsed(maxIterationsUsed),
      _timeWindowsReused(timeWindowsReused),
      _primaryDataIDs(std::move(dataIDs)),
      _onBoundViolation(std::move(onBoundViolation)),
      _forceInitialRelaxation(forceInitialRelaxation),
      _reducedTimeGrid(reducedTimeGrid),
      _qrV(filter),
      _filter(filter),
      _singularityLimit(singularityLimit),
      _infostringstream(std::ostringstream::ate)
{
  PRECICE_CHECK((_initialRelaxation > 0.0) && (_initialRelaxation <= 1.0),
                "Initial relaxation factor for QN acceleration has to "
                "be larger than zero and smaller or equal than one. "
                "Current initial relaxation is {}",
                _initialRelaxation);
  PRECICE_CHECK(_maxIterationsUsed > 0,
                "Maximum number of iterations used in the quasi-Newton acceleration "
                "scheme has to be larger than zero. "
                "Current maximum reused iterations is {}",
                _maxIterationsUsed);
  PRECICE_CHECK(_timeWindowsReused >= 0,
                "Number of previous time windows to be reused for "
                "quasi-Newton acceleration has to be larger than or equal to zero. "
                "Current number of time windows reused is {}",
                _timeWindowsReused);
}

/** ---------------------------------------------------------------------------------------------
 *         initialize()
 *
 * @brief: Initializes all the needed variables and data
 *  ---------------------------------------------------------------------------------------------
 */
void BaseQNAcceleration::initialize(
    const DataMap &cplData)
{
  PRECICE_TRACE(cplData.size());

  for (const DataMap::value_type &pair : cplData) {
    PRECICE_ASSERT(pair.second->getSize() == pair.second->getPreviousIterationSize(), "current and previousIteration have to be initialized and of identical size.",
                   pair.second->getSize(), pair.second->getPreviousIterationSize());
  }

  PRECICE_WARN_IF(
      std::any_of(cplData.cbegin(), cplData.cend(), [](const auto &p) { return p.second->hasGradient(); }),
      "Gradient data, which is required by at least one of the configured data mappings, is not yet compatible with quasi-Newton acceleration. This combination might lead to numerical issues. "
      "Consider switching to a different acceleration scheme or a different data mapping scheme.");

  checkDataIDs(cplData);
  // store all data IDs in vector
  _dataIDs.clear();
  for (const DataMap::value_type &pair : cplData) {
    _dataIDs.push_back(pair.first);
  }

  _matrixCols.clear();
  _matrixCols.push_front(0);
  _firstIteration  = true;
  _firstTimeWindow = true;
}

void BaseQNAcceleration::checkBound(Eigen::VectorXd &data, DataMap &cplData, const std::vector<DataID> &dataIDs, std::string onBoundViolation, Eigen::VectorXd &xUpdate)
{
  Eigen::Index offset            = 0;
  bool         violationDetected = false;
  for (auto id : dataIDs) {
    if (violationDetected)
      break;

    Eigen::Index size          = cplData.at(id)->values().size();
    int          dataDimension = cplData.at(id)->getDimensions();
    auto         lowerBound    = cplData.at(id)->getLowerBound();
    auto         upperBound    = cplData.at(id)->getUpperBound();

    for (int j = 0; j < dataDimension; j++) {
      for (Eigen::Index i = j; i < size; i += dataDimension) {
        if (data[i + offset] < lowerBound[j] || data[i + offset] > upperBound[j]) {
          violationDetected = true;
          break;
        }
      }
    }
    offset += size;
  }
  if (violationDetected) {
    if (onBoundViolation == "discard") {
      PRECICE_WARN("The coupling data has violated its bound after the Quasi-Newton step. The current step will be discarded.");
      data -= xUpdate;
    } else if (onBoundViolation == "clamp") {
      PRECICE_WARN("The coupling data has violated its bound after the Quasi-Newton step. The values will be clamped to their bounds.");

      offset = 0;
      for (auto id : dataIDs) {
        Eigen::Index size          = cplData.at(id)->values().size();
        int          dataDimension = cplData.at(id)->getDimensions();
        auto         lowerBound    = cplData.at(id)->getLowerBound();
        auto         upperBound    = cplData.at(id)->getUpperBound();
        for (int j = 0; j < dataDimension; j++) {
          for (Eigen::Index i = j; i < size; i += dataDimension) {
            data[i + offset] = std::clamp(data[i + offset], lowerBound[j], upperBound[j]);
          }
        }
        offset += size;
      }
    } else if (onBoundViolation == "scale") {
      PRECICE_WARN(
          "The coupling data has violated its bound after the Quasi-Newton step. The step length will be scaled to avoid the bound violation.");
      offset           = 0;
      double scaleStep = 0.0;
      for (auto id : dataIDs) {
        Eigen::Index size          = cplData.at(id)->values().size();
        int          dataDimension = cplData.at(id)->getDimensions();
        auto         lowerBound    = cplData.at(id)->getLowerBound();
        auto         upperBound    = cplData.at(id)->getUpperBound();

        for (int j = 0; j < dataDimension; j++) {
          for (Eigen::Index i = j; i < size; i += dataDimension) {
            if (xUpdate[i + offset] > 0 && data[i + offset] > upperBound[j])
              scaleStep = std::max(scaleStep, (data[i + offset] - upperBound[j]) / xUpdate[i + offset]);
            else if (xUpdate[i + offset] < 0 && data[i + offset] < lowerBound[j])
              scaleStep = std::max(scaleStep, (data[i + offset] - lowerBound[j]) / xUpdate[i + offset]);
          }
        }
        offset += size;
      }
      data -= xUpdate * scaleStep;
    }
  }
}
/** ---------------------------------------------------------------------------------------------
 *         updateDifferenceMatrices()
 *
 * @brief: computes the current residual and stores it, computes the differences and
 *         updates the difference matrices F and C.
 *  ---------------------------------------------------------------------------------------------
 */
void BaseQNAcceleration::updateDifferenceMatrices(
    const DataMap &cplData)
{
  PRECICE_TRACE();

  // Compute current residual: vertex-data - oldData
  _primaryResiduals = _primaryValues;
  _primaryResiduals -= _oldPrimaryValues;
  _residuals = _values;
  _residuals -= _oldValues;

  PRECICE_WARN_IF(math::equals(utils::IntraComm::l2norm(_primaryResiduals), 0.0),
                  "The coupling residual equals almost zero. There is maybe something wrong in your adapter. "
                  "Maybe you always write the same data or you call advance without "
                  "providing new data first or you do not use available read data. "
                  "Or you just converge much further than actually necessary.");

  // if (_firstIteration && (_firstTimeWindow || (_matrixCols.size() < 2))) {
  if (_firstIteration && (_firstTimeWindow || _forceInitialRelaxation)) {

    _aitkenFactor = _initialRelaxation;

    // do nothing: constant relaxation
  } else {
    PRECICE_DEBUG("   Update Difference Matrices");
    if (not _firstIteration) {
      // Update matrices V, W with newest information

      PRECICE_ASSERT(_matrixV.cols() == _matrixW.cols(), _matrixV.cols(), _matrixW.cols());
      PRECICE_ASSERT(getLSSystemCols() <= _maxIterationsUsed, getLSSystemCols(), _maxIterationsUsed);

      PRECICE_WARN_IF(
          2 * getLSSystemCols() >= getLSSystemRows(),
          "The number of columns in the least squares system exceeded half the number of unknowns at the interface. "
          "The system will probably become bad or ill-conditioned and the quasi-Newton acceleration may not "
          "converge. Maybe the number of allowed columns (\"max-used-iterations\") should be limited.");

      Eigen::VectorXd deltaR = _primaryResiduals;
      deltaR -= _oldPrimaryResiduals;

      Eigen::VectorXd deltaXTilde = _values;
      deltaXTilde -= _oldXTilde;

      double residualMagnitude = utils::IntraComm::l2norm(deltaR);

      if (not math::equals(utils::IntraComm::l2norm(_primaryValues), 0.0)) {
        residualMagnitude /= utils::IntraComm::l2norm(_primaryValues);
      }
      PRECICE_WARN_IF(
          math::equals(residualMagnitude, 0.0),
          "Adding a vector with a two-norm of {} to the quasi-Newton V matrix, which will lead to "
          "ill-conditioning. A filter might delete the column again. Still, this could mean that you are "
          "converging too tightly, that you reached steady-state, or that you are giving by mistake identical "
          "data to preCICE in two consecutive iterations.",
          residualMagnitude);

      bool columnLimitReached = getLSSystemCols() == _maxIterationsUsed;
      bool overdetermined     = getLSSystemCols() <= getLSSystemRows();
      if (not columnLimitReached && overdetermined) {

        utils::appendFront(_matrixV, deltaR);
        utils::appendFront(_matrixW, deltaXTilde);

        // insert column deltaR = _primaryResiduals - _oldPrimaryResiduals at pos. 0 (front) into the
        // QR decomposition and update decomposition

        // apply scaling here
        _preconditioner->apply(deltaR);
        _qrV.pushFront(deltaR);

        _matrixCols.front()++;
      } else {
        utils::shiftSetFirst(_matrixV, deltaR);
        utils::shiftSetFirst(_matrixW, deltaXTilde);

        // inserts column deltaR at pos. 0 to the QR decomposition and deletes the last column
        // the QR decomposition of V is updated
        _preconditioner->apply(deltaR);
        _qrV.pushFront(deltaR);
        _qrV.popBack();

        _matrixCols.front()++;
        _matrixCols.back()--;
        if (_matrixCols.back() == 0) {
          _matrixCols.pop_back();
        }
        _nbDropCols++;
      }
    }
    _oldPrimaryResiduals = _primaryResiduals; // Store residuals
    _oldPrimaryXTilde    = _primaryValues;    // Store x_tilde
    _oldXTilde           = _values;           // Store coupling x_tilde
  }
}

/** ---------------------------------------------------------------------------------------------
 *         performAcceleration()
 *
 * @brief: performs one iteration of the quasi Newton acceleration.
 *  ---------------------------------------------------------------------------------------------
 */
void BaseQNAcceleration::performAcceleration(
    DataMap &cplData,
    double   windowStart,
    double   windowEnd)
{
  PRECICE_TRACE(_primaryDataIDs.size(), cplData.size());

  profiling::Event e("cpl.computeQuasiNewtonUpdate", profiling::Synchronize);

  // We can only initialize here since we need to know the time grid already.
  if (_firstTimeWindow and _firstIteration) {
    initializeVectorsAndPreconditioner(cplData, windowStart);
  }

  PRECICE_ASSERT(_oldPrimaryResiduals.size() == _oldPrimaryXTilde.size(), _oldPrimaryResiduals.size(), _oldPrimaryXTilde.size());
  PRECICE_ASSERT(_primaryValues.size() == _oldPrimaryXTilde.size(), _primaryValues.size(), _oldPrimaryXTilde.size());
  PRECICE_ASSERT(_oldPrimaryValues.size() == _oldPrimaryXTilde.size(), _oldPrimaryValues.size(), _oldPrimaryXTilde.size());
  PRECICE_ASSERT(_primaryResiduals.size() == _oldPrimaryXTilde.size(), _primaryResiduals.size(), _oldPrimaryXTilde.size());

  // Needs to be called in the first iteration of each time window.
  if (_firstIteration) {
    _timeGrids.value().moveTimeGridToNewWindow(cplData);
    _primaryTimeGrids.value().moveTimeGridToNewWindow(cplData);
  }

  /// Sample all the data to the corresponding time grid in _timeGrids and concatenate everything into a long vector
  /// timeGrids are stored using std::optional, thus the .value() to get the actual object
  concatenateCouplingData(_values, _oldValues, cplData, _dataIDs, _timeGrids.value(), windowStart);
  concatenateCouplingData(_primaryValues, _oldPrimaryValues, cplData, _primaryDataIDs, _primaryTimeGrids.value(), windowStart);

  /** update the difference matrices V,W  includes:
   * scaling of values
   * computation of residuals
   * appending the difference matrices
   */
  updateDifferenceMatrices(cplData);

  if (_firstIteration && (_firstTimeWindow || _forceInitialRelaxation)) {
    PRECICE_DEBUG("   Performing underrelaxation");
    _oldPrimaryXTilde    = _primaryValues;    // Store x tilde of primary data
    _oldXTilde           = _values;           // Store x tilde of primary and secondary data
    _oldPrimaryResiduals = _primaryResiduals; // Store current residual of primary data

    applyRelaxation(_initialRelaxation, cplData, windowStart);
    its++;
    _firstIteration = false;
    return;
  }

  PRECICE_DEBUG("   Performing quasi-Newton Step");

  // If the previous time window converged within one single iteration, nothing was added
  // to the LS system matrices and they need to be restored from the backup at time T-2
  if (not _firstTimeWindow && (getLSSystemCols() < 1) && (_timeWindowsReused == 0) && not _forceInitialRelaxation) {
    PRECICE_DEBUG("   Last time window converged after one iteration. Need to restore the matrices from backup.");

    _matrixCols = _matrixColsBackup;
    _matrixV    = _matrixVBackup;
    _matrixW    = _matrixWBackup;

    // re-computation of QR decomposition from _matrixV = _matrixVBackup
    // this occurs very rarely, to be precise, it occurs only if the coupling terminates
    // after the first iteration and the matrix data from time window t-2 has to be used
    _preconditioner->apply(_matrixV);
    _qrV.reset(_matrixV, getLSSystemRows());
    _preconditioner->revert(_matrixV);
    _resetLS = true; // need to recompute _Wtil, Q, R (only for IMVJ efficient update)
  }

  /**
   *  === update and apply preconditioner ===
   *
   * The preconditioner is only applied to the matrix V and the columns that are inserted into the
   * QR-decomposition of V.
   */
  _preconditioner->update(false, _primaryValues, _primaryResiduals);

  // apply scaling to V, V' := P * V (only needed to reset the QR-dec of V)
  _preconditioner->apply(_matrixV);

  if (_preconditioner->requireNewQR()) {
    if (not(_filter == Acceleration::QR2FILTER || _filter == Acceleration::QR3FILTER)) { // for QR2 and QR3 filter, there is no need to do this twice
      _qrV.reset(_matrixV, getLSSystemRows());
    }
    if (_filter == Acceleration::QR3FILTER) { // QR3 filter needs to recompute QR3 filter
      _qrV.requireQR3Fallback();
    }
    _preconditioner->newQRfulfilled();
  }

  if (_firstIteration) {
    _nbDelCols  = 0;
    _nbDropCols = 0;
  }

  // apply the configured filter to the LS system
  profiling::Event applyingFilter("ApplyFilter");
  applyFilter();
  applyingFilter.stop();

  // revert scaling of V, in computeQNUpdate all data objects are unscaled.
  _preconditioner->revert(_matrixV);

  /**
   * compute quasi-Newton update
   * PRECONDITION: All objects are unscaled, except the matrices within the QR-dec of V.
   *               Thus, the pseudo inverse needs to be reverted before using it.
   */
  Eigen::VectorXd xUpdate = Eigen::VectorXd::Zero(_values.size());
  computeQNUpdate(xUpdate);

  // Apply the quasi-Newton update
  _values += xUpdate;

  // Check for bound violations
  if (_onBoundViolation != "ignore")
    checkBound(_values, cplData, _dataIDs, _onBoundViolation, xUpdate);

  // pending deletion: delete old V, W matrices if timeWindowsReused = 0
  // those were only needed for the first iteration (instead of underrelax.)
  if (_firstIteration && _timeWindowsReused == 0 && not _forceInitialRelaxation) {
    // save current matrix data in case the coupling for the next time window will terminate
    // after the first iteration (no new data, i.e., V = W = 0)
    if (getLSSystemCols() > 0) {
      _matrixColsBackup = _matrixCols;
      _matrixVBackup    = _matrixV;
      _matrixWBackup    = _matrixW;
    }
    // if no time windows reused, the matrix data needs to be cleared as it was only needed for the
    // QN-step in the first iteration (idea: rather perform QN-step with information from last converged
    // time window instead of doing a underrelaxation)
    if (not _firstTimeWindow) {
      _matrixV.resize(0, 0);
      _matrixW.resize(0, 0);
      _matrixCols.clear();
      _matrixCols.push_front(0); // vital after clear()
      _qrV.reset();
      // set the number of global rows in the QRFactorization.
      _qrV.setGlobalRows(getPrimaryLSSystemRows());
      _resetLS = true; // need to recompute _Wtil, Q, R (only for IMVJ efficient update)
    }
  }

  PRECICE_CHECK(
      !std::isnan(utils::IntraComm::l2norm(xUpdate)),
      "The quasi-Newton update contains NaN values. This means that the quasi-Newton acceleration failed to converge. "
      "When writing your own adapter this could indicate that you give wrong information to preCICE, such as identical "
      "data in succeeding iterations. Or you do not properly save and reload checkpoints. "
      "If you give the correct data this could also mean that the coupled problem is too hard to solve. Try to use a QR "
      "filter or increase its threshold (larger epsilon).");

  // updates the waveform and values in coupling data by splitting the primary and secondary data back into the individual data objects
  updateCouplingData(cplData, windowStart);

  // number of iterations (usually equals number of columns in LS-system)
  its++;
  _firstIteration = false;
}

void BaseQNAcceleration::applyFilter()
{
  PRECICE_TRACE(_filter);

  if (_filter == Acceleration::NOFILTER) {
    // do nothing
  } else {
    // do: filtering of least-squares system to maintain good conditioning
    std::vector<int> delIndices(0);
    _qrV.applyFilter(_singularityLimit, delIndices, _matrixV);
    // start with largest index (as V,W matrices are shrunk and shifted

    for (int i = delIndices.size() - 1; i >= 0; i--) {

      removeMatrixColumn(delIndices[i]);

      PRECICE_DEBUG(" Filter: removing column with index {} in iteration {} of time window: {}", delIndices[i], its, tWindows);
    }
    PRECICE_ASSERT(_matrixV.cols() == _qrV.cols(), _matrixV.cols(), _qrV.cols());
  }
}

void BaseQNAcceleration::updateCouplingData(
    const DataMap &cplData, double windowStart)
{

  PRECICE_TRACE();
  // offset to keep track of the position in xUpdate
  Eigen::Index offset = 0;

  for (int id : _dataIDs) {

    auto  &couplingData = *cplData.at(id);
    size_t dataSize     = couplingData.getSize();

    Eigen::VectorXd timeGrid = _timeGrids->getTimeGridAfter(id, windowStart);
    couplingData.waveform().trimAfter(windowStart);
    for (int i = 0; i < timeGrid.size(); i++) {

      Eigen::VectorXd temp = Eigen::VectorXd::Zero(dataSize);
      for (size_t j = 0; j < dataSize; j++) {
        temp(j) = _values(offset + j);
      }
      offset += dataSize;

      couplingData.setSampleAtTime(timeGrid(i), time::Sample(couplingData.getDimensions(), temp));
    }
  }
}

/** ---------------------------------------------------------------------------------------------
 *         iterationsConverged()
 *
 * @brief: Is called when the convergence criterion for the coupling is fulfilled and finalizes
 *         the quasi Newton acceleration. Stores new differences in F and C, clears or
 *         updates F and C according to the number of reused time windows
 *  ---------------------------------------------------------------------------------------------
 */
void BaseQNAcceleration::iterationsConverged(
    const DataMap &cplData, double windowStart)
{
  PRECICE_TRACE();

  if (utils::IntraComm::isPrimary() || !utils::IntraComm::isParallel())
    _infostringstream << "# time window " << tWindows << " converged #\n iterations: " << its
                      << "\n used cols: " << getLSSystemCols() << "\n del cols: " << _nbDelCols << '\n';

  its = 0;
  tWindows++;

  // the most recent differences for the V, W matrices have not been added so far
  // this has to be done in iterations converged, as PP won't be called any more if
  // convergence was achieved

  // If, in the first time window, we converge already in the first iteration, we have not yet initialized. Then, we need to do it here.
  if (_firstTimeWindow and _firstIteration) {
    initializeVectorsAndPreconditioner(cplData, windowStart);
  }

  // Needs to be called in the first iteration of each time window.
  if (_firstIteration) {
    _timeGrids->moveTimeGridToNewWindow(cplData);
    _primaryTimeGrids->moveTimeGridToNewWindow(cplData);
  }
  /// Sample all the data to the corresponding time grid in _timeGrids and concatenate everything into a long vector
  /// timeGrids are stored using std::optional, thus the .value() to get the actual object
  concatenateCouplingData(_values, _oldValues, cplData, _dataIDs, _timeGrids.value(), windowStart);
  concatenateCouplingData(_primaryValues, _oldPrimaryValues, cplData, _primaryDataIDs, _primaryTimeGrids.value(), windowStart);
  updateDifferenceMatrices(cplData);

  if (not _matrixCols.empty() && _matrixCols.front() == 0) { // Did only one iteration
    _matrixCols.pop_front();
  }

#ifndef NDEBUG
  std::ostringstream stream;
  stream << "Matrix column counters: ";
  for (int cols : _matrixCols) {
    stream << cols << ", ";
  }
  PRECICE_DEBUG(stream.str());
#endif // Debug

  // doing specialized stuff for the corresponding acceleration scheme after
  // convergence of iteration i.e.:
  // - analogously to the V,W matrices, remove columns from matrices for secondary data
  // - save the old Jacobian matrix
  specializedIterationsConverged(cplData);

  // if we already have convergence in the first iteration of the first time window
  // we need to do underrelaxation in the first iteration of the second time window
  // so "_firstTimeWindow" is slightly misused, but still the best way to understand
  // the concept
  if (not _firstIteration)
    _firstTimeWindow = false;

  // update preconditioner depending on residuals or values (must be after specialized iterations converged --> IMVJ)
  _preconditioner->update(true, _primaryValues, _primaryResiduals);

  if (_timeWindowsReused == 0) {
    if (_forceInitialRelaxation) {
      _matrixV.resize(0, 0);
      _matrixW.resize(0, 0);
      _qrV.reset();
      // set the number of global rows in the QRFactorization.
      _qrV.setGlobalRows(getPrimaryLSSystemRows());
      _matrixCols.clear(); // _matrixCols.push_front() at the end of the method.
    } else {
      /**
       * pending deletion (after first iteration of next time window
       * Using the matrices from the old time window for the first iteration
       * is better than doing underrelaxation as first iteration of every time window
       */
    }
  } else if (static_cast<int>(_matrixCols.size()) > _timeWindowsReused) {
    int toRemove = _matrixCols.back();
    _nbDropCols += toRemove;
    PRECICE_ASSERT(toRemove > 0, toRemove);
    PRECICE_DEBUG("Removing {} cols from least-squares system with {} cols", toRemove, getLSSystemCols());
    PRECICE_ASSERT(_matrixV.cols() == _matrixW.cols(), _matrixV.cols(), _matrixW.cols());
    PRECICE_ASSERT(getLSSystemCols() > toRemove, getLSSystemCols(), toRemove);

    // remove columns
    for (int i = 0; i < toRemove; i++) {
      utils::removeColumnFromMatrix(_matrixV, _matrixV.cols() - 1);
      utils::removeColumnFromMatrix(_matrixW, _matrixW.cols() - 1);
      // also remove the corresponding columns from the dynamic QR-descomposition of _matrixV
      _qrV.popBack();
    }
    _matrixCols.pop_back();
  }

  _matrixCols.push_front(0);
  _firstIteration = true;
}

/** ---------------------------------------------------------------------------------------------
 *         removeMatrixColumn()
 *
 * @brief: removes a column from the least squares system, i. e., from the matrices F and C
 *  ---------------------------------------------------------------------------------------------
 */
void BaseQNAcceleration::removeMatrixColumn(
    int columnIndex)
{
  PRECICE_TRACE(columnIndex, _matrixV.cols());

  _nbDelCols++;

  PRECICE_ASSERT(_matrixV.cols() > 1);
  utils::removeColumnFromMatrix(_matrixV, columnIndex);
  utils::removeColumnFromMatrix(_matrixW, columnIndex);

  // Reduce column count
  std::deque<int>::iterator iter = _matrixCols.begin();
  int                       cols = 0;
  while (iter != _matrixCols.end()) {
    cols += *iter;
    if (cols > columnIndex) {
      PRECICE_ASSERT(*iter > 0);
      *iter -= 1;
      if (*iter == 0) {
        _matrixCols.erase(iter);
      }
      break;
    }
    iter++;
  }
}

void BaseQNAcceleration::exportState(
    io::TXTWriter &writer)
{
}

void BaseQNAcceleration::importState(
    io::TXTReader &reader)
{
}

int BaseQNAcceleration::getDeletedColumns() const
{
  return _nbDelCols;
}

int BaseQNAcceleration::getDroppedColumns() const
{
  return _nbDropCols;
}

int BaseQNAcceleration::getLSSystemCols() const
{
  int cols = 0;
  for (int col : _matrixCols) {
    cols += col;
  }
  if (_hasNodesOnInterface) {
    PRECICE_ASSERT(cols == _matrixV.cols(), cols, _matrixV.cols(), _matrixCols, _qrV.cols());
    PRECICE_ASSERT(cols == _matrixW.cols(), cols, _matrixW.cols());
  }

  return cols;
}

int BaseQNAcceleration::getMaxUsedIterations() const
{
  return _maxIterationsUsed;
}

int BaseQNAcceleration::getMaxUsedTimeWindows() const
{
  return _timeWindowsReused;
}

int BaseQNAcceleration::getLSSystemRows() const
{
  if (utils::IntraComm::isParallel()) {
    return _dimOffsets.back();
  }
  return _residuals.size();
}

int BaseQNAcceleration::getPrimaryLSSystemRows() const
{
  if (utils::IntraComm::isParallel()) {
    return _dimOffsetsPrimary.back();
  }
  return _primaryResiduals.size();
}

void BaseQNAcceleration::writeInfo(
    const std::string &s, bool allProcs)
{
  if (not utils::IntraComm::isParallel()) {
    // serial acceleration mode
    _infostringstream << s;

    // parallel acceleration
  } else {
    if (not allProcs) {
      if (utils::IntraComm::isPrimary())
        _infostringstream << s;
    } else {
      _infostringstream << s;
    }
  }
  _infostringstream << std::flush;
}

void BaseQNAcceleration::concatenateCouplingData(Eigen::VectorXd &data, Eigen::VectorXd &oldData, const DataMap &cplData, std::vector<int> dataIDs, precice::time::TimeGrids timeGrids, double windowStart) const
{
  Eigen::Index offset = 0;
  for (int id : dataIDs) {
    Eigen::Index    dataSize = cplData.at(id)->getSize();
    Eigen::VectorXd timeGrid = timeGrids.getTimeGridAfter(id, windowStart);

    for (int i = 0; i < timeGrid.size(); i++) {

      auto current = cplData.at(id)->waveform().sample(timeGrid(i));
      auto old     = cplData.at(id)->getPreviousValuesAtTime(timeGrid(i));

      PRECICE_ASSERT(data.size() >= offset + dataSize, "the values were not initialized correctly");
      PRECICE_ASSERT(oldData.size() >= offset + dataSize, "the values were not initialized correctly");

      for (Eigen::Index i = 0; i < dataSize; i++) {
        data(i + offset)    = current(i);
        oldData(i + offset) = old(i);
      }
      offset += dataSize;
    }
  }
}

void BaseQNAcceleration::initializeVectorsAndPreconditioner(const DataMap &cplData, double windowStart)
{
  // Saves the time grid of each waveform in the data field to be used in the QN method
  _timeGrids.emplace(cplData, _dataIDs, false);
  _primaryTimeGrids.emplace(cplData, _primaryDataIDs, _reducedTimeGrid);

  // Helper function
  auto addTimeSliceSize = [&](size_t sum, int id, precice::time::TimeGrids timeGrids) { return sum + timeGrids.getTimeGridAfter(id, windowStart).size() * cplData.at(id)->getSize(); };

  // Size of primary data
  const size_t primaryDataSize = std::accumulate(_primaryDataIDs.begin(), _primaryDataIDs.end(), (size_t) 0, [&](size_t sum, int id) { return addTimeSliceSize(sum, id, _primaryTimeGrids.value()); });

  // Size of values
  const size_t dataSize = std::accumulate(_dataIDs.begin(), _dataIDs.end(), (size_t) 0, [&](size_t sum, int id) { return addTimeSliceSize(sum, id, _timeGrids.value()); });

  _values              = Eigen::VectorXd::Zero(dataSize);
  _oldValues           = Eigen::VectorXd::Zero(dataSize);
  _oldXTilde           = Eigen::VectorXd::Zero(dataSize);
  _residuals           = Eigen::VectorXd::Zero(dataSize);
  _primaryValues       = Eigen::VectorXd::Zero(primaryDataSize);
  _oldPrimaryValues    = Eigen::VectorXd::Zero(primaryDataSize);
  _oldPrimaryXTilde    = Eigen::VectorXd::Zero(primaryDataSize);
  _primaryResiduals    = Eigen::VectorXd::Zero(primaryDataSize);
  _oldPrimaryResiduals = Eigen::VectorXd::Zero(primaryDataSize);

  // Also clear the matrices, since the initialization phase will be called more than once when using remeshing
  _matrixW.resize(0, 0);
  _matrixV.resize(0, 0);

  /**
   *  make dimensions public to all procs,
   *  last entry _dimOffsets[IntraComm::getSize()] holds the global dimension, global,n
   */
  std::stringstream ss;
  if (utils::IntraComm::isParallel()) {
    PRECICE_ASSERT(utils::IntraComm::getCommunication() != nullptr);
    PRECICE_ASSERT(utils::IntraComm::getCommunication()->isConnected());

    if (primaryDataSize <= 0) {
      _hasNodesOnInterface = false;
    }

    /** provide vertex offset information for all processors
     *  mesh->getVertexOffsets() provides an array that stores the number of mesh vertices on each processor
     *  This information needs to be gathered for all meshes. To get the number of respective unknowns of a specific processor
     *  we need to multiply the number of vertices with the dimensionality of the vector-valued data for each coupling data.
     */
    _dimOffsets.resize(utils::IntraComm::getSize() + 1);
    _dimOffsetsPrimary.resize(utils::IntraComm::getSize() + 1);
    _dimOffsets[0]        = 0;
    _dimOffsetsPrimary[0] = 0;
    for (size_t i = 0; i < _dimOffsets.size() - 1; i++) {
      int accumulatedNumberOfUnknowns        = 0;
      int accumulatedNumberOfPrimaryUnknowns = 0;

      for (auto &elem : _dataIDs) {
        const auto &offsets = cplData.at(elem)->getVertexOffsets();

        accumulatedNumberOfUnknowns += offsets[i] * cplData.at(elem)->getDimensions() * _timeGrids.value().getTimeGridAfter(elem, windowStart).size();

        if (utils::contained(elem, _primaryDataIDs)) {
          accumulatedNumberOfPrimaryUnknowns += offsets[i] * cplData.at(elem)->getDimensions() * _primaryTimeGrids.value().getTimeGridAfter(elem, windowStart).size();
        }
      }
      _dimOffsets[i + 1]        = accumulatedNumberOfUnknowns;
      _dimOffsetsPrimary[i + 1] = accumulatedNumberOfPrimaryUnknowns;
    }

    PRECICE_DEBUG("Number of unknowns at the interface (global): {}", _dimOffsets.back());
    if (utils::IntraComm::isPrimary()) {
      _infostringstream << fmt::format("\n--------\n DOFs (global): {}\n offsets: {}\n", _dimOffsets.back(), _dimOffsets);
    }

    // test that the computed number of unknown per proc equals the number of primaryDataSize actually present on that proc
    const size_t unknowns = _dimOffsets[utils::IntraComm::getRank() + 1] - _dimOffsets[utils::IntraComm::getRank()];
    PRECICE_ASSERT(dataSize == unknowns, dataSize, unknowns);
    const size_t primaryUnknowns = _dimOffsetsPrimary[utils::IntraComm::getRank() + 1] - _dimOffsetsPrimary[utils::IntraComm::getRank()];
    PRECICE_ASSERT(primaryDataSize == primaryUnknowns, primaryDataSize, primaryUnknowns);
  } else {
    _infostringstream << fmt::format("\n--------\n DOFs (global): {}\n", dataSize);
  }

  // set the number of global rows in the QRFactorization.
  _qrV.reset();
  _qrV.setGlobalRows(getPrimaryLSSystemRows());

  std::vector<size_t> subVectorSizes; // needed for preconditioner
  std::transform(_primaryDataIDs.cbegin(), _primaryDataIDs.cend(), std::back_inserter(subVectorSizes), [&cplData, windowStart, this](const auto &d) { return _primaryTimeGrids.value().getTimeGridAfter(d, windowStart).size() * cplData.at(d)->getSize(); });
  _preconditioner->initialize(subVectorSizes);

  specializedInitializeVectorsAndPreconditioner(cplData);
}

} // namespace acceleration
} // namespace precice
