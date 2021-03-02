#include "acceleration/BaseQNAcceleration.hpp"
#include <Eigen/Core>
#include <cmath>
#include <memory>
#include "acceleration/impl/Preconditioner.hpp"
#include "acceleration/impl/QRFactorization.hpp"
#include "com/Communication.hpp"
#include "com/SharedPointer.hpp"
#include "cplscheme/CouplingData.hpp"
#include "logging/LogMacros.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "utils/EigenHelperFunctions.hpp"
#include "utils/Event.hpp"
#include "utils/Helpers.hpp"
#include "utils/MasterSlave.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace io {
class TXTReader;
class TXTWriter;
} // namespace io

extern bool syncMode;
namespace acceleration {

/* ----------------------------------------------------------------------------
 *     Constructor
 * ----------------------------------------------------------------------------
 */
BaseQNAcceleration::BaseQNAcceleration(
    double                  initialRelaxation,
    bool                    forceInitialRelaxation,
    int                     maxIterationsUsed,
    int                     timestepsReused,
    int                     filter,
    double                  singularityLimit,
    std::vector<int>        dataIDs,
    impl::PtrPreconditioner preconditioner)
    : _preconditioner(preconditioner),
      _initialRelaxation(initialRelaxation),
      _maxIterationsUsed(maxIterationsUsed),
      _timestepsReused(timestepsReused),
      _dataIDs(dataIDs),
      _forceInitialRelaxation(forceInitialRelaxation),
      _qrV(filter),
      _filter(filter),
      _singularityLimit(singularityLimit),
      _infostringstream(std::ostringstream::ate)
{
  PRECICE_CHECK((_initialRelaxation > 0.0) && (_initialRelaxation <= 1.0),
                "Initial relaxation factor for QN acceleration has to "
                    << "be larger than zero and smaller or equal than one. Current initial relaxation is: " << _initialRelaxation);
  PRECICE_CHECK(_maxIterationsUsed > 0,
                "Maximum number of iterations used in the quasi-Newton acceleration "
                    << "scheme has to be larger than zero. Current maximum reused iterations is: " << _maxIterationsUsed);
  PRECICE_CHECK(_timestepsReused >= 0,
                "Number of previous time windows to be reused for quasi-Newton acceleration has to be larger than or equal to zero. "
                    << "Current number of time windows reused is " << _timestepsReused);
}

/** ---------------------------------------------------------------------------------------------
 *         initialize()
 *
 * @brief: Initializes all the needed variables and data
 *  ---------------------------------------------------------------------------------------------
 */
void BaseQNAcceleration::initialize(
    DataMap &cplData)
{
  PRECICE_TRACE(cplData.size());
  checkDataIDs(cplData);

  /*
  std::stringstream sss;
  sss<<"debugOutput-rank-"<<utils::MasterSlave::getRank();
  _debugOut.open(sss.str(), std::ios_base::out);
  _debugOut << std::setprecision(16);

  Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", " << ", ";");

  _debugOut<<"initialization:\n";
  for (int id : _dataIDs) {
      const auto& values = *cplData[id]->values;
      const auto& oldValues = cplData[id]->oldValues.col(0);

      _debugOut<<"id: "<<id<<" dim: "<<cplData[id]->dimension<<"     values: "<<values.format(CommaInitFmt)<<'\n';
      _debugOut<<"id: "<<id<<" dim: "<<cplData[id]->dimension<<" old values: "<<oldValues.format(CommaInitFmt)<<'\n';
    }
  _debugOut<<"\n";
  */

  size_t              entries = 0;
  std::vector<size_t> subVectorSizes; //needed for preconditioner

  for (auto &elem : _dataIDs) {
    entries += cplData[elem]->values().size();
    subVectorSizes.push_back(cplData[elem]->values().size());
  }

  _matrixCols.push_front(0);
  _firstIteration = true;
  _firstTimeStep  = true;

  PRECICE_ASSERT(_oldXTilde.size() == 0);
  PRECICE_ASSERT(_oldResiduals.size() == 0);
  _oldXTilde    = Eigen::VectorXd::Zero(entries);
  _oldResiduals = Eigen::VectorXd::Zero(entries);
  _residuals    = Eigen::VectorXd::Zero(entries);
  _values       = Eigen::VectorXd::Zero(entries);
  _oldValues    = Eigen::VectorXd::Zero(entries);

  /**
   *  make dimensions public to all procs,
   *  last entry _dimOffsets[MasterSlave::getSize()] holds the global dimension, global,n
   */
  std::stringstream ss;
  if (utils::MasterSlave::isMaster() || utils::MasterSlave::isSlave()) {
    PRECICE_ASSERT(utils::MasterSlave::_communication.get() != NULL);
    PRECICE_ASSERT(utils::MasterSlave::_communication->isConnected());

    if (entries <= 0) {
      _hasNodesOnInterface = false;
    }

    /** provide vertex offset information for all processors
     *  mesh->getVertexOffsets() provides an array that stores the number of mesh vertices on each processor
     *  This information needs to be gathered for all meshes. To get the number of respective unknowns of a specific processor
     *  we need to multiply the number of vertices with the dimensionality of the vector-valued data for each coupling data.
     */
    _dimOffsets.resize(utils::MasterSlave::getSize() + 1);
    _dimOffsets[0] = 0;
    //for (auto & elem : _dataIDs) {
    //	std::cout<<" Offsets:(vertex) \n"<<cplData[elem]->mesh->getVertexOffsets()<<'\n';
    //}
    for (size_t i = 0; i < _dimOffsets.size() - 1; i++) {
      int accumulatedNumberOfUnknowns = 0;
      for (auto &elem : _dataIDs) {
        auto &offsets = cplData[elem]->mesh->getVertexOffsets();
        accumulatedNumberOfUnknowns += offsets[i] * cplData[elem]->getDimensions();
      }
      _dimOffsets[i + 1] = accumulatedNumberOfUnknowns;
    }
    PRECICE_DEBUG("Number of unknowns at the interface (global): " << _dimOffsets.back());
    if (utils::MasterSlave::isMaster()) {
      _infostringstream << "\n--------\n DOFs (global): " << _dimOffsets.back() << "\n offsets: " << _dimOffsets << '\n';
    }

    // test that the computed number of unknown per proc equals the number of entries actually present on that proc
    size_t unknowns = _dimOffsets[utils::MasterSlave::getRank() + 1] - _dimOffsets[utils::MasterSlave::getRank()];
    PRECICE_ASSERT(entries == unknowns, entries, unknowns);
  } else {
    _infostringstream << "\n--------\n DOFs (global): " << entries << '\n';
  }

  // set the number of global rows in the QRFactorization. This is essential for the correctness in master-slave mode!
  _qrV.setGlobalRows(getLSSystemRows());

  // Fetch secondary data IDs, to be relaxed with same coefficients from IQN-ILS
  for (DataMap::value_type &pair : cplData) {
    if (not utils::contained(pair.first, _dataIDs)) {
      _secondaryDataIDs.push_back(pair.first);
      int secondaryEntries            = pair.second->values().size();
      _secondaryResiduals[pair.first] = Eigen::VectorXd::Zero(secondaryEntries);
    }
  }

  // Append old value columns, if not done outside of acceleration already
  for (DataMap::value_type &pair : cplData) {
    int cols = pair.second->oldValues.cols();
    if (cols < 1) { // Add only, if not already done
      //PRECICE_ASSERT(pair.second->values().size() > 0, pair.first);
      utils::append(pair.second->oldValues, (Eigen::VectorXd) Eigen::VectorXd::Zero(pair.second->values().size()));
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
    DataMap &cplData)
{
  PRECICE_TRACE();

  // Compute current residual: vertex-data - oldData
  _residuals = _values;
  _residuals -= _oldValues;

  if (math::equals(utils::MasterSlave::l2norm(_residuals), 0.0)) {
    PRECICE_WARN("The coupling residual equals almost zero. There is maybe something wrong in your adapter. "
                 "Maybe you always write the same data or you call advance without "
                 "providing new data first or you do not use available read data. "
                 "Or you just converge much further than actually necessary.");
  }

  //if (_firstIteration && (_firstTimeStep || (_matrixCols.size() < 2))) {
  if (_firstIteration && (_firstTimeStep || _forceInitialRelaxation)) {
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
            << "The system will probably become bad or ill-conditioned and the quasi-Newton acceleration may not "
            << "converge. Maybe the number of allowed columns (\"max-used-iterations\") should be limited.");

      Eigen::VectorXd deltaR = _residuals;
      deltaR -= _oldResiduals;

      Eigen::VectorXd deltaXTilde = _values;
      deltaXTilde -= _oldXTilde;

      PRECICE_CHECK(not math::equals(utils::MasterSlave::l2norm(deltaR), 0.0), "Attempting to add a zero vector to the quasi-Newton V matrix. This means that the residual "
                                                                               "in two consecutive iterations is identical. There is probably something wrong in your adapter. "
                                                                               "Maybe you always write the same (or only incremented) data or you call advance without "
                                                                               "providing  new data first.");

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
    DataMap &cplData)
{
  PRECICE_TRACE(_dataIDs.size(), cplData.size());

  utils::Event e("cpl.computeQuasiNewtonUpdate", precice::syncMode);

  PRECICE_ASSERT(_oldResiduals.size() == _oldXTilde.size(), _oldResiduals.size(), _oldXTilde.size());
  PRECICE_ASSERT(_values.size() == _oldXTilde.size(), _values.size(), _oldXTilde.size());
  PRECICE_ASSERT(_oldValues.size() == _oldXTilde.size(), _oldValues.size(), _oldXTilde.size());
  PRECICE_ASSERT(_residuals.size() == _oldXTilde.size(), _residuals.size(), _oldXTilde.size());

  /*
  Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", " << ", ";");
  _debugOut<<"iteration: "<<its<<" tStep: "<<tSteps<<"   cplData entry:\n";
  for (int id : _dataIDs) {
      const auto& values = *cplData[id]->values;
      const auto& oldValues = cplData[id]->oldValues.col(0);

      _debugOut<<"id: "<<id<<"     values: "<<values.format(CommaInitFmt)<<'\n';
      _debugOut<<"id: "<<id<<" old values: "<<oldValues.format(CommaInitFmt)<<'\n';
    }
  _debugOut<<"\n";
  */

  // assume data structures associated with the LS system can be updated easily.

  // scale data values (and secondary data values)
  concatenateCouplingData(cplData);

  /** update the difference matrices V,W  includes:
   * scaling of values
   * computation of residuals
   * appending the difference matrices
   */
  updateDifferenceMatrices(cplData);

  if (_firstIteration && (_firstTimeStep || _forceInitialRelaxation)) {
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

    // If the previous time step converged within one single iteration, nothing was added
    // to the LS system matrices and they need to be restored from the backup at time T-2
    if (not _firstTimeStep && (getLSSystemCols() < 1) && (_timestepsReused == 0) && not _forceInitialRelaxation) {
      PRECICE_DEBUG("   Last time step converged after one iteration. Need to restore the matrices from backup.");

      _matrixCols = _matrixColsBackup;
      _matrixV    = _matrixVBackup;
      _matrixW    = _matrixWBackup;

      // re-computation of QR decomposition from _matrixV = _matrixVBackup
      // this occurs very rarely, to be precise, it occurs only if the coupling terminates
      // after the first iteration and the matrix data from time step t-2 has to be used
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
      if (not(_filter == Acceleration::QR2FILTER)) { //for QR2 filter, there is no need to do this twice
        _qrV.reset(_matrixV, getLSSystemRows());
      }
      _preconditioner->newQRfulfilled();
    }

    if (_firstIteration) {
      _nbDelCols  = 0;
      _nbDropCols = 0;
    }

    // apply the configured filter to the LS system
    applyFilter();

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

    // pending deletion: delete old V, W matrices if timestepsReused = 0
    // those were only needed for the first iteration (instead of underrelax.)
    if (_firstIteration && _timestepsReused == 0 && not _forceInitialRelaxation) {
      // save current matrix data in case the coupling for the next time step will terminate
      // after the first iteration (no new data, i.e., V = W = 0)
      if (getLSSystemCols() > 0) {
        _matrixColsBackup = _matrixCols;
        _matrixVBackup    = _matrixV;
        _matrixWBackup    = _matrixW;
      }
      // if no time steps reused, the matrix data needs to be cleared as it was only needed for the
      // QN-step in the first iteration (idea: rather perform QN-step with information from last converged
      // time step instead of doing a underrelaxation)
      if (not _firstTimeStep) {
        _matrixV.resize(0, 0);
        _matrixW.resize(0, 0);
        _matrixCols.clear();
        _matrixCols.push_front(0); // vital after clear()
        _qrV.reset();
        // set the number of global rows in the QRFactorization. This is essential for the correctness in master-slave mode!
        _qrV.setGlobalRows(getLSSystemRows());
        _resetLS = true; // need to recompute _Wtil, Q, R (only for IMVJ efficient update)
      }
    }

    if (std::isnan(utils::MasterSlave::l2norm(xUpdate))) {
      PRECICE_ERROR("The quasi-Newton update contains NaN values. This means that the quasi-Newton acceleration failed to converge. "
                    "When writing your own adapter this could indicate that you give wrong information to preCICE, such as identical "
                    "data in succeeding iterations. Or you do not properly save and reload checkpoints. "
                    "If you give the correct data this could also mean that the coupled problem is too hard to solve. Try to use a QR "
                    "filter or increase its threshold (larger epsilon).");
    }
  }

  splitCouplingData(cplData);

  /*
  _debugOut<<"finished update: \n";
  for (int id : _dataIDs) {
      const auto& values = *cplData[id]->values;
      const auto& oldValues = cplData[id]->oldValues.col(0);

      _debugOut<<"id: "<<id<<"norm: "<<values.norm()<<"     values: "<<values.format(CommaInitFmt)<<'\n';
      _debugOut<<"id: "<<id<<"norm: "<<oldValues.norm()<<" old values: "<<oldValues.format(CommaInitFmt)<<'\n';
    }
  _debugOut<<"\n";
  */

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
    // start with largest index (as V,W matrices are shrinked and shifted

    for (int i = delIndices.size() - 1; i >= 0; i--) {

      removeMatrixColumn(delIndices[i]);

      PRECICE_DEBUG(" Filter: removing column with index " << delIndices[i] << " in iteration " << its << " of time step: " << tSteps);
    }
    PRECICE_ASSERT(_matrixV.cols() == _qrV.cols(), _matrixV.cols(), _qrV.cols());
  }
}

void BaseQNAcceleration::concatenateCouplingData(
    DataMap &cplData)
{
  PRECICE_TRACE();

  int offset = 0;
  for (int id : _dataIDs) {
    int         size      = cplData[id]->values().size();
    auto &      values    = cplData[id]->values();
    const auto &oldValues = cplData[id]->oldValues.col(0);
    for (int i = 0; i < size; i++) {
      _values(i + offset)    = values(i);
      _oldValues(i + offset) = oldValues(i);
    }
    offset += size;
  }
}

void BaseQNAcceleration::splitCouplingData(
    DataMap &cplData)
{
  PRECICE_TRACE();

  int offset = 0;
  for (int id : _dataIDs) {
    int   size       = cplData[id]->values().size();
    auto &valuesPart = cplData[id]->values();
    //Eigen::VectorXd& oldValuesPart = cplData[id]->oldValues.col(0);
    cplData[id]->oldValues.col(0) = _oldValues.segment(offset, size); /// @todo: check if this is correct
    for (int i = 0; i < size; i++) {
      valuesPart(i) = _values(i + offset);
      //oldValuesPart(i) = _oldValues(i + offset);
    }
    offset += size;
  }
}

/** ---------------------------------------------------------------------------------------------
 *         iterationsConverged()
 *
 * @brief: Is called when the convergence criterion for the coupling is fullfilled and finalizes
 *         the quasi Newton acceleration. Stores new differences in F and C, clears or
 *         updates F and C according to the number of reused time steps
 *  ---------------------------------------------------------------------------------------------
 */
void BaseQNAcceleration::iterationsConverged(
    DataMap &cplData)
{
  PRECICE_TRACE();

  if (utils::MasterSlave::isMaster() || (not utils::MasterSlave::isMaster() && not utils::MasterSlave::isSlave()))
    _infostringstream << "# time step " << tSteps << " converged #\n iterations: " << its
                      << "\n used cols: " << getLSSystemCols() << "\n del cols: " << _nbDelCols << '\n';

  its = 0;
  tSteps++;

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

  // if we already have convergence in the first iteration of the first timestep
  // we need to do underrelax in the first iteration of the second timesteps
  // so "_firstTimeStep" is slightly misused, but still the best way to understand
  // the concept
  if (not _firstIteration)
    _firstTimeStep = false;

  // update preconditioner depending on residuals or values (must be after specialized iterations converged --> IMVJ)
  _preconditioner->update(true, _values, _residuals);

  if (_timestepsReused == 0) {
    if (_forceInitialRelaxation) {
      _matrixV.resize(0, 0);
      _matrixW.resize(0, 0);
      _qrV.reset();
      // set the number of global rows in the QRFactorization. This is essential for the correctness in master-slave mode!
      _qrV.setGlobalRows(getLSSystemRows());
      _matrixCols.clear(); // _matrixCols.push_front() at the end of the method.
    } else {
      /**
       * pending deletion (after first iteration of next time step
       * Using the matrices from the old time step for the first iteration
       * is better than doing underrelaxation as first iteration of every time step
       */
    }
  } else if ((int) _matrixCols.size() > _timestepsReused) {
    int toRemove = _matrixCols.back();
    _nbDropCols += toRemove;
    PRECICE_ASSERT(toRemove > 0, toRemove);
    PRECICE_DEBUG("Removing " << toRemove << " cols from least-squares system with " << getLSSystemCols() << " cols");
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
  if (utils::MasterSlave::isMaster() || utils::MasterSlave::isSlave()) {
    return _dimOffsets.back();
  }
  return _residuals.size();
}

void BaseQNAcceleration::writeInfo(
    std::string s, bool allProcs)
{
  if (not utils::MasterSlave::isMaster() && not utils::MasterSlave::isSlave()) {
    // serial acceleration mode
    _infostringstream << s;

    // parallel acceleration, master-slave mode
  } else {
    if (not allProcs) {
      if (utils::MasterSlave::isMaster())
        _infostringstream << s;
    } else {
      _infostringstream << s;
    }
  }
  _infostringstream << std::flush;
}
} // namespace acceleration
} // namespace precice
