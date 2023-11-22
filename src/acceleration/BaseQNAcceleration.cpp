#include "acceleration/BaseQNAcceleration.hpp"
#include <Eigen/Core>
#include <boost/range/adaptor/map.hpp>
#include <cmath>
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
    impl::PtrPreconditioner preconditioner)
    : _preconditioner(std::move(preconditioner)),
      _initialRelaxation(initialRelaxation),
      _maxIterationsUsed(maxIterationsUsed),
      _timeWindowsReused(timeWindowsReused),
      _dataIDs(std::move(dataIDs)),
      _forceInitialRelaxation(forceInitialRelaxation),
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

  if (std::any_of(cplData.cbegin(), cplData.cend(), [](const auto &p) { return p.second->hasGradient(); })) {
    PRECICE_WARN("Gradient data, which is required by at least one of the configured data mappings, is not yet compatible with quasi-Newton acceleration. This combination might lead to numerical issues. "
                 "Consider switching to a different acceleration scheme or a different data mapping scheme.");
  }

  checkDataIDs(cplData);

  for (const auto &data : cplData | boost::adaptors::map_values) {
    if (data->exchangeSubsteps()) {
      PRECICE_ERROR("Quasi-Newton acceleration does not yet support using data from all substeps. Please set substeps=\"false\" in the exchange tag of data \"{}\".", data->getDataName());
    }
  }

  size_t              entries = 0;
  std::vector<size_t> subVectorSizes; // needed for preconditioner

  for (auto &elem : _dataIDs) {
    entries += cplData.at(elem)->getSize();
    subVectorSizes.push_back(cplData.at(elem)->getSize());
  }

  _matrixCols.push_front(0);
  _firstIteration  = true;
  _firstTimeWindow = true;

  PRECICE_ASSERT(_oldXTilde.size() == 0);
  PRECICE_ASSERT(_oldResiduals.size() == 0);
  _oldXTilde    = Eigen::VectorXd::Zero(entries);
  _oldResiduals = Eigen::VectorXd::Zero(entries);
  _residuals    = Eigen::VectorXd::Zero(entries);
  _values       = Eigen::VectorXd::Zero(entries);
  _oldValues    = Eigen::VectorXd::Zero(entries);

  /**
   *  make dimensions public to all procs,
   *  last entry _dimOffsets[IntraComm::getSize()] holds the global dimension, global,n
   */
  std::stringstream ss;
  if (utils::IntraComm::isParallel()) {
    PRECICE_ASSERT(utils::IntraComm::getCommunication() != nullptr);
    PRECICE_ASSERT(utils::IntraComm::getCommunication()->isConnected());

    if (entries <= 0) {
      _hasNodesOnInterface = false;
    }

    /** provide vertex offset information for all processors
     *  mesh->getVertexOffsets() provides an array that stores the number of mesh vertices on each processor
     *  This information needs to be gathered for all meshes. To get the number of respective unknowns of a specific processor
     *  we need to multiply the number of vertices with the dimensionality of the vector-valued data for each coupling data.
     */
    _dimOffsets.resize(utils::IntraComm::getSize() + 1);
    _dimOffsets[0] = 0;
    for (size_t i = 0; i < _dimOffsets.size() - 1; i++) {
      int accumulatedNumberOfUnknowns = 0;
      for (auto &elem : _dataIDs) {
        const auto &offsets = cplData.at(elem)->getVertexOffsets();
        accumulatedNumberOfUnknowns += offsets[i] * cplData.at(elem)->getDimensions();
      }
      _dimOffsets[i + 1] = accumulatedNumberOfUnknowns;
    }
    PRECICE_DEBUG("Number of unknowns at the interface (global): {}", _dimOffsets.back());
    if (utils::IntraComm::isPrimary()) {
      _infostringstream << fmt::format("\n--------\n DOFs (global): {}\n offsets: {}\n", _dimOffsets.back(), _dimOffsets);
    }

    // test that the computed number of unknown per proc equals the number of entries actually present on that proc
    size_t unknowns = _dimOffsets[utils::IntraComm::getRank() + 1] - _dimOffsets[utils::IntraComm::getRank()];
    PRECICE_ASSERT(entries == unknowns, entries, unknowns);
  } else {
    _infostringstream << fmt::format("\n--------\n DOFs (global): {}\n", entries);
  }

  // set the number of global rows in the QRFactorization.
  _qrV.setGlobalRows(getLSSystemRows());

  // Fetch secondary data IDs, to be relaxed with same coefficients from IQN-ILS
  for (const DataMap::value_type &pair : cplData) {
    if (not utils::contained(pair.first, _dataIDs)) {
      _secondaryDataIDs.push_back(pair.first);
      int secondaryEntries            = pair.second->getSize();
      _secondaryResiduals[pair.first] = Eigen::VectorXd::Zero(secondaryEntries);
    }
  }

  _preconditioner->initialize(subVectorSizes);
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
  _residuals = _values;
  _residuals -= _oldValues;

  if (math::equals(utils::IntraComm::l2norm(_residuals), 0.0)) {
    PRECICE_WARN("The coupling residual equals almost zero. There is maybe something wrong in your adapter. "
                 "Maybe you always write the same data or you call advance without "
                 "providing new data first or you do not use available read data. "
                 "Or you just converge much further than actually necessary.");
  }

  // if (_firstIteration && (_firstTimeWindow || (_matrixCols.size() < 2))) {
  if (_firstIteration && (_firstTimeWindow || _forceInitialRelaxation)) {
    // do nothing: constant relaxation
  } else {
    PRECICE_DEBUG("   Update Difference Matrices");
    if (not _firstIteration) {
      // Update matrices V, W with newest information

      PRECICE_ASSERT(_matrixV.cols() == _matrixW.cols(), _matrixV.cols(), _matrixW.cols());
      PRECICE_ASSERT(getLSSystemCols() <= _maxIterationsUsed, getLSSystemCols(), _maxIterationsUsed);

      if (2 * getLSSystemCols() >= getLSSystemRows())
        PRECICE_WARN(
            "The number of columns in the least squares system exceeded half the number of unknowns at the interface. "
            "The system will probably become bad or ill-conditioned and the quasi-Newton acceleration may not "
            "converge. Maybe the number of allowed columns (\"max-used-iterations\") should be limited.");

      Eigen::VectorXd deltaR = _residuals;
      deltaR -= _oldResiduals;

      Eigen::VectorXd deltaXTilde = _values;
      deltaXTilde -= _oldXTilde;

      double residualMagnitude = utils::IntraComm::l2norm(deltaR);

      if (not math::equals(utils::IntraComm::l2norm(_values), 0.0)) {
        residualMagnitude /= utils::IntraComm::l2norm(_values);
      }
      if (math::equals(residualMagnitude, 0.0)) {
        PRECICE_WARN("Adding a vector with a two-norm of {} to the quasi-Newton V matrix, which will lead to "
                     "ill-conditioning. A filter might delete the column again. Still, this could mean that you are " 
                     "converging too tightly, that you reached steady-state, or that you are giving by mistake identical "
                     "data to preCICE in two consecutive iterations.", residualMagnitude);
      }

      bool columnLimitReached = getLSSystemCols() == _maxIterationsUsed;
      bool overdetermined     = getLSSystemCols() <= getLSSystemRows();
      if (not columnLimitReached && overdetermined) {

        utils::appendFront(_matrixV, deltaR);
        utils::appendFront(_matrixW, deltaXTilde);

        // insert column deltaR = _residuals - _oldResiduals at pos. 0 (front) into the
        // QR decomposition and update decomposition

        //apply scaling here
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
    _oldResiduals = _residuals; // Store residuals
    _oldXTilde    = _values;    // Store x_tilde
  }
}

/** ---------------------------------------------------------------------------------------------
 *         performAcceleration()
 *
 * @brief: performs one iteration of the quasi Newton acceleration.
 *  ---------------------------------------------------------------------------------------------
 */
void BaseQNAcceleration::performAcceleration(
    const DataMap &cplData)
{
  PRECICE_TRACE(_dataIDs.size(), cplData.size());

  profiling::Event e("cpl.computeQuasiNewtonUpdate", profiling::Synchronize);

  PRECICE_ASSERT(_oldResiduals.size() == _oldXTilde.size(), _oldResiduals.size(), _oldXTilde.size());
  PRECICE_ASSERT(_values.size() == _oldXTilde.size(), _values.size(), _oldXTilde.size());
  PRECICE_ASSERT(_oldValues.size() == _oldXTilde.size(), _oldValues.size(), _oldXTilde.size());
  PRECICE_ASSERT(_residuals.size() == _oldXTilde.size(), _residuals.size(), _oldXTilde.size());

  // assume data structures associated with the LS system can be updated easily.

  // scale data values (and secondary data values)
  concatenateCouplingData(cplData);

  /** update the difference matrices V,W  includes:
   * scaling of values
   * computation of residuals
   * appending the difference matrices
   */
  updateDifferenceMatrices(cplData);

  if (_firstIteration && (_firstTimeWindow || _forceInitialRelaxation)) {
    PRECICE_DEBUG("   Performing underrelaxation");
    _oldXTilde    = _values;    // Store x tilde
    _oldResiduals = _residuals; // Store current residual

    // Perform constant relaxation
    // with residual: x_new = x_old + omega * res
    _residuals *= _initialRelaxation;
    _residuals += _oldValues;
    _values = _residuals;

    computeUnderrelaxationSecondaryData(cplData);
  } else {
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

    _preconditioner->update(false, _values, _residuals);
    // apply scaling to V, V' := P * V (only needed to reset the QR-dec of V)
    _preconditioner->apply(_matrixV);

    if (_preconditioner->requireNewQR()) {
      if (not(_filter == Acceleration::QR2FILTER)) { // for QR2 filter, there is no need to do this twice
        _qrV.reset(_matrixV, getLSSystemRows());
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
    Eigen::VectorXd xUpdate = Eigen::VectorXd::Zero(_residuals.size());
    computeQNUpdate(cplData, xUpdate);

    /**
     * apply quasiNewton update
     */
    _values = _oldValues + xUpdate + _residuals; // = x^k + delta_x + r^k - q^k

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
        _qrV.setGlobalRows(getLSSystemRows());
        _resetLS = true; // need to recompute _Wtil, Q, R (only for IMVJ efficient update)
      }
    }

    if (std::isnan(utils::IntraComm::l2norm(xUpdate))) {
      PRECICE_ERROR("The quasi-Newton update contains NaN values. This means that the quasi-Newton acceleration failed to converge. "
                    "When writing your own adapter this could indicate that you give wrong information to preCICE, such as identical "
                    "data in succeeding iterations. Or you do not properly save and reload checkpoints. "
                    "If you give the correct data this could also mean that the coupled problem is too hard to solve. Try to use a QR "
                    "filter or increase its threshold (larger epsilon).");
    }
  }

  splitCouplingData(cplData);
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

void BaseQNAcceleration::concatenateCouplingData(
    const DataMap &cplData)
{
  PRECICE_TRACE();

  int offset = 0;
  for (int id : _dataIDs) {
    int         size      = cplData.at(id)->getSize();
    auto &      values    = cplData.at(id)->values();
    const auto &oldValues = cplData.at(id)->previousIteration();
    for (int i = 0; i < size; i++) {
      _values(i + offset)    = values(i);
      _oldValues(i + offset) = oldValues(i);
    }
    offset += size;
  }
}

void BaseQNAcceleration::splitCouplingData(
    const DataMap &cplData)
{
  PRECICE_TRACE();

  int offset = 0;
  for (int id : _dataIDs) {
    int   size       = cplData.at(id)->getSize();
    auto &valuesPart = cplData.at(id)->values();
    for (int i = 0; i < size; i++) {
      valuesPart(i) = _values(i + offset);
    }
    offset += size;
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
    const DataMap &cplData)
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
  concatenateCouplingData(cplData);
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
  _preconditioner->update(true, _values, _residuals);

  if (_timeWindowsReused == 0) {
    if (_forceInitialRelaxation) {
      _matrixV.resize(0, 0);
      _matrixW.resize(0, 0);
      _qrV.reset();
      // set the number of global rows in the QRFactorization.
      _qrV.setGlobalRows(getLSSystemRows());
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

int BaseQNAcceleration::getLSSystemRows()
{
  if (utils::IntraComm::isParallel()) {
    return _dimOffsets.back();
  }
  return _residuals.size();
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
} // namespace acceleration
} // namespace precice
