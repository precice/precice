#include "BaseQNPostProcessing.hpp"
#include <sstream>
#include "QRFactorization.hpp"
#include "com/Communication.hpp"
#include "cplscheme/CouplingData.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "utils/EigenHelperFunctions.hpp"
#include "utils/MasterSlave.hpp"
#include "utils/Helpers.hpp"
#include "utils/Event.hpp"

namespace precice
{
extern bool syncMode;
namespace cplscheme
{
namespace impl
{

/* ----------------------------------------------------------------------------
 *     Constructor
 * ----------------------------------------------------------------------------
 */
BaseQNPostProcessing::BaseQNPostProcessing(
    double            initialRelaxation,
    bool              forceInitialRelaxation,
    int               maxIterationsUsed,
    int               timestepsReused,
    int               filter,
    double            singularityLimit,
    std::vector<int>  dataIDs,
    PtrPreconditioner preconditioner)
  :   _preconditioner(preconditioner),
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
  CHECK((_initialRelaxation > 0.0) && (_initialRelaxation <= 1.0),
        "Initial relaxation factor for QN post-processing has to "
            << "be larger than zero and smaller or equal than one!");
  CHECK(_maxIterationsUsed > 0,
        "Maximal iterations used for QN post-processing has to be larger than zero!");
  CHECK(_timestepsReused >= 0,
        "Number of old timesteps to be reused for QN post-processing has to be >= 0!");
}

/** ---------------------------------------------------------------------------------------------
 *         initialize()
 *
 * @brief: Initializes all the needed variables and data
 *  ---------------------------------------------------------------------------------------------
 */
void BaseQNPostProcessing::initialize(
    DataMap &cplData)
{
  TRACE(cplData.size());

  /*
  std::stringstream sss;
  sss<<"debugOutput-rank-"<<utils::MasterSlave::_rank;
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
    CHECK(utils::contained(elem, cplData),
          "Data with ID " << elem << " is not contained in data given at initialization!");
    entries += cplData[elem]->values->size();
    subVectorSizes.push_back(cplData[elem]->values->size());
  }

  _matrixCols.push_front(0);
  _firstIteration = true;
  _firstTimeStep  = true;

  assertion(_oldXTilde.size() == 0);
  assertion(_oldResiduals.size() == 0);
  _oldXTilde    = Eigen::VectorXd::Zero(entries);
  _oldResiduals = Eigen::VectorXd::Zero(entries);
  _residuals    = Eigen::VectorXd::Zero(entries);
  _values       = Eigen::VectorXd::Zero(entries);
  _oldValues    = Eigen::VectorXd::Zero(entries);

  // if design specifiaction not initialized yet
  if (not(_designSpecification.size() > 0)) {
    _designSpecification = Eigen::VectorXd::Zero(_residuals.size());
  }
  /**
   *  make dimensions public to all procs,
   *  last entry _dimOffsets[MasterSlave::_size] holds the global dimension, global,n
   */
  std::stringstream ss;
  if (utils::MasterSlave::_masterMode || utils::MasterSlave::_slaveMode) {
    assertion(utils::MasterSlave::_communication.get() != NULL);
    assertion(utils::MasterSlave::_communication->isConnected());

    if (entries <= 0) {
      _hasNodesOnInterface = false;
    }

    /** provide vertex offset information for all processors
     *  mesh->getVertexOffsets() provides an array that stores the number of mesh vertices on each processor
     *  This information needs to be gathered for all meshes. To get the number of respective unknowns of a specific processor
     *  we need to multiply the number of vertices with the dimensionality of the vector-valued data for each coupling data.
     */
    _dimOffsets.resize(utils::MasterSlave::_size + 1);
    _dimOffsets[0] = 0;
    //for (auto & elem : _dataIDs) {
    //	std::cout<<" Offsets:(vertex) \n"<<cplData[elem]->mesh->getVertexOffsets()<<'\n';
    //}
    for (size_t i = 0; i < _dimOffsets.size() - 1; i++) {
      int accumulatedNumberOfUnknowns = 0;
      for (auto &elem : _dataIDs) {
        auto &offsets = cplData[elem]->mesh->getVertexOffsets();
        accumulatedNumberOfUnknowns += offsets[i] * cplData[elem]->dimension;
      }
      _dimOffsets[i + 1] = accumulatedNumberOfUnknowns;
    }
    DEBUG("Number of unknowns at the interface (global): " << _dimOffsets.back());
    if (utils::MasterSlave::_masterMode) {
      _infostringstream << "\n--------\n DOFs (global): " << _dimOffsets.back() << "\n offsets: " << _dimOffsets << '\n';
    }

    // test that the computed number of unknown per proc equals the number of entries actually present on that proc
    size_t unknowns = _dimOffsets[utils::MasterSlave::_rank + 1] - _dimOffsets[utils::MasterSlave::_rank];
    assertion(entries == unknowns, entries, unknowns);
  } else {
    _infostringstream << "\n--------\n DOFs (global): " << entries << '\n';
  }

  // set the number of global rows in the QRFactorization. This is essential for the correctness in master-slave mode!
  _qrV.setGlobalRows(getLSSystemRows());

  // Fetch secondary data IDs, to be relaxed with same coefficients from IQN-ILS
  for (DataMap::value_type &pair : cplData) {
    if (not utils::contained(pair.first, _dataIDs)) {
      _secondaryDataIDs.push_back(pair.first);
      int secondaryEntries            = pair.second->values->size();
      _secondaryResiduals[pair.first] = Eigen::VectorXd::Zero(secondaryEntries);
    }
  }

  // Append old value columns, if not done outside of post-processing already
  for (DataMap::value_type &pair : cplData) {
    int cols = pair.second->oldValues.cols();
    if (cols < 1) { // Add only, if not already done
      //assertion(pair.second->values->size() > 0, pair.first);
      utils::append(pair.second->oldValues, (Eigen::VectorXd) Eigen::VectorXd::Zero(pair.second->values->size()));
    }
  }

  _preconditioner->initialize(subVectorSizes);
}

/** ---------------------------------------------------------------------------------------------
 *         setDesignSpecification()
 *
 * @brief: sets a design specification for the fine model optimization problem
 *         i.e., x_star = argmin_x || f(x) - q ||
 *  ---------------------------------------------------------------------------------------------
 */
void BaseQNPostProcessing::setDesignSpecification(
    Eigen::VectorXd &q)
{
  TRACE();
  assertion(q.size() == _residuals.size(), q.size(), _residuals.size());
  _designSpecification = q;
}

/** ---------------------------------------------------------------------------------------------
 *         getDesignSpecification()
 *
 * @brief: Returns the design specification corresponding to the given coupling data.
 *         This information is needed for convergence measurements in the coupling scheme.
 *  ---------------------------------------------------------------------------------------------
 */
std::map<int, Eigen::VectorXd> BaseQNPostProcessing::getDesignSpecification(
    DataMap &cplData)
{
  TRACE();
  std::map<int, Eigen::VectorXd> designSpecifications;
  int                            off = 0;
  for (int id : _dataIDs) {
    int             size = cplData[id]->values->size();
    Eigen::VectorXd q    = Eigen::VectorXd::Zero(size);
    for (int i = 0; i < size; i++) {
      q(i) = _designSpecification(i + off);
    }
    off += size;
    std::map<int, Eigen::VectorXd>::value_type pair = std::make_pair(id, q);
    designSpecifications.insert(pair);
  }
  return designSpecifications;
}

/** ---------------------------------------------------------------------------------------------
 *         updateDifferenceMatrices()
 *
 * @brief: computes the current coarse and fine model residual, computes the differences and
 *         updates the difference matrices F and C. Also stores the residuals
 *  ---------------------------------------------------------------------------------------------
 */
void BaseQNPostProcessing::updateDifferenceMatrices(
    DataMap &cplData)
{
  TRACE();
  
  // Compute current residual: vertex-data - oldData
  _residuals = _values;
  _residuals -= _oldValues;

  //if (_firstIteration && (_firstTimeStep || (_matrixCols.size() < 2))) {
  if (_firstIteration && (_firstTimeStep || _forceInitialRelaxation)) {
    // do nothing: constant relaxation
  } else {
    DEBUG("   Update Difference Matrices");
    if (not _firstIteration) {
      // Update matrices V, W with newest information

      assertion(_matrixV.cols() == _matrixW.cols(), _matrixV.cols(), _matrixW.cols());
      assertion(getLSSystemCols() <= _maxIterationsUsed, getLSSystemCols(), _maxIterationsUsed);

      if (2 * getLSSystemCols() >= getLSSystemRows())
        WARN(
            "The number of columns in the least squares system exceeded half the number of unknowns at the interface. The system will probably become bad or ill-conditioned and the quasi-Newton post processing may not converge. Maybe the number of allowed columns (maxIterationsUsed) should be limited.");

      Eigen::VectorXd deltaR = _residuals;
      deltaR -= _oldResiduals;

      Eigen::VectorXd deltaXTilde = _values;
      deltaXTilde -= _oldXTilde;

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
      }
    }
    _oldResiduals = _residuals; // Store residuals
    _oldXTilde    = _values;    // Store x_tilde
  }
}

/** ---------------------------------------------------------------------------------------------
 *         performPostProcessing()
 *
 * @brief: performs one iteration of the quasi Newton post processing. It steers the execution
 *         of fine and coarse model evaluations and also calls the coarse model optimization.
 *  ---------------------------------------------------------------------------------------------
 */
void BaseQNPostProcessing::performPostProcessing(
    DataMap &cplData)
{
  TRACE(_dataIDs.size(), cplData.size());
  
  utils::Event e("cpl.computeQuasiNewtonUpdate", precice::syncMode);

  assertion(_oldResiduals.size() == _oldXTilde.size(), _oldResiduals.size(), _oldXTilde.size());
  assertion(_values.size() == _oldXTilde.size(), _values.size(), _oldXTilde.size());
  assertion(_oldValues.size() == _oldXTilde.size(), _oldValues.size(), _oldXTilde.size());
  assertion(_residuals.size() == _oldXTilde.size(), _residuals.size(), _oldXTilde.size());

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
    DEBUG("   Performing underrelaxation");
    _oldXTilde    = _values;    // Store x tilde
    _oldResiduals = _residuals; // Store current residual

    // Perform constant relaxation
    // with residual: x_new = x_old + omega * (res-q)
    _residuals *= _initialRelaxation;
    _residuals -= (_designSpecification * _initialRelaxation);
    _residuals += _oldValues;
    _values = _residuals;

    computeUnderrelaxationSecondaryData(cplData);
  } else {
    DEBUG("   Performing quasi-Newton Step");

    // If the previous time step converged within one single iteration, nothing was added
    // to the LS system matrices and they need to be restored from the backup at time T-2
    if (not _firstTimeStep && (getLSSystemCols() < 1) && (_timestepsReused == 0) && not _forceInitialRelaxation) {
      DEBUG("   Last time step converged after one iteration. Need to restore the matrices from backup.");

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

    // subtract design specification from residuals, i.e., we want to minimize argmin_x|| r(x) - q ||
    assertion(_residuals.size() == _designSpecification.size(), _residuals.size(), _designSpecification.size());
    _residuals -= _designSpecification;

    /**
     *  === update and apply preconditioner ===
     *
     * The preconditioner is only applied to the matrix V and the columns that are inserted into the
     * QR-decomposition of V.
     * Note: here, the _residuals are H(x)- x - q, i.e., residual of the fixed-point iteration
     *       minus the design specification of the optimization problem (!= null if MM is used)
     */

    _preconditioner->update(false, _values, _residuals);
    // apply scaling to V, V' := P * V (only needed to reset the QR-dec of V)
    _preconditioner->apply(_matrixV);

    if (_preconditioner->requireNewQR()) {
      if (not(_filter == PostProcessing::QR2FILTER)) { //for QR2 filter, there is no need to do this twice
        _qrV.reset(_matrixV, getLSSystemRows());
      }
      _preconditioner->newQRfulfilled();
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

    /// todo maybe add design specification. Though, residuals are overwritten in the next iteration this would be a clearer and nicer code

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
      ERROR("The coupling iteration in time step " << tSteps << " failed to converge and NaN values occurred throughout the coupling process. ");
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

void BaseQNPostProcessing::applyFilter()
{
  TRACE(_filter);
  
  if (_filter == PostProcessing::NOFILTER) {
    // do nothing
  } else {
    // do: filtering of least-squares system to maintain good conditioning
    std::vector<int> delIndices(0);
    _qrV.applyFilter(_singularityLimit, delIndices, _matrixV);
    // start with largest index (as V,W matrices are shrinked and shifted
    for (int i = delIndices.size() - 1; i >= 0; i--) {

      removeMatrixColumn(delIndices[i]);

      DEBUG(" Filter: removing column with index " << delIndices[i] << " in iteration " << its << " of time step: " << tSteps);
    }
    assertion(_matrixV.cols() == _qrV.cols(), _matrixV.cols(), _qrV.cols());
  }
}

void BaseQNPostProcessing::concatenateCouplingData(
    DataMap &cplData)
{
  TRACE();

  int offset = 0;
  for (int id : _dataIDs) {
    int         size      = cplData[id]->values->size();
    auto &      values    = *cplData[id]->values;
    const auto &oldValues = cplData[id]->oldValues.col(0);
    for (int i = 0; i < size; i++) {
      _values(i + offset)    = values(i);
      _oldValues(i + offset) = oldValues(i);
    }
    offset += size;
  }
}

void BaseQNPostProcessing::splitCouplingData(
    DataMap &cplData)
{
  TRACE();

  int offset = 0;
  for (int id : _dataIDs) {
    int   size       = cplData[id]->values->size();
    auto &valuesPart = *(cplData[id]->values);
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
 *         the quasi Newton post processing. Stores new differences in F and C, clears or
 *         updates F and C according to the number of reused time steps
 *  ---------------------------------------------------------------------------------------------
 */
void BaseQNPostProcessing::iterationsConverged(
    DataMap &cplData)
{
  TRACE();
  
  if (utils::MasterSlave::_masterMode || (not utils::MasterSlave::_masterMode && not utils::MasterSlave::_slaveMode))
    _infostringstream << "# time step " << tSteps << " converged #\n iterations: " << its
                      << "\n used cols: " << getLSSystemCols() << "\n del cols: " << _nbDelCols << '\n';

  its = 0;
  tSteps++;
  _nbDelCols = 0;

  // the most recent differences for the V, W matrices have not been added so far
  // this has to be done in iterations converged, as PP won't be called any more if
  // convergence was achieved
  concatenateCouplingData(cplData);
  updateDifferenceMatrices(cplData);

  // subtract design specification from residuals, i.e., we want to minimize argmin_x|| r(x) - q ||
  assertion(_residuals.size() == _designSpecification.size(), _residuals.size(), _designSpecification.size());
  _residuals -= _designSpecification;

  if (_matrixCols.front() == 0) { // Did only one iteration
    _matrixCols.pop_front();
  }

#ifndef NDEBUG
  std::ostringstream stream;
  stream << "Matrix column counters: ";
  for (int cols : _matrixCols) {
    stream << cols << ", ";
  }
  DEBUG(stream.str());
#endif // Debug

  // doing specialized stuff for the corresponding post processing scheme after
  // convergence of iteration i.e.:
  // - analogously to the V,W matrices, remove columns from matrices for secondary data
  // - save the old Jacobian matrix
  specializedIterationsConverged(cplData);

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
    assertion(toRemove > 0, toRemove);
    DEBUG("Removing " << toRemove << " cols from least-squares system with " << getLSSystemCols() << " cols");
    assertion(_matrixV.cols() == _matrixW.cols(), _matrixV.cols(), _matrixW.cols());
    assertion(getLSSystemCols() > toRemove, getLSSystemCols(), toRemove);

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
void BaseQNPostProcessing::removeMatrixColumn(
    int columnIndex)
{
  TRACE(columnIndex, _matrixV.cols());

  // debugging information, can be removed
  _nbDelCols++;

  assertion(_matrixV.cols() > 1);
  utils::removeColumnFromMatrix(_matrixV, columnIndex);
  utils::removeColumnFromMatrix(_matrixW, columnIndex);

  // Reduce column count
  std::deque<int>::iterator iter = _matrixCols.begin();
  int                       cols = 0;
  while (iter != _matrixCols.end()) {
    cols += *iter;
    if (cols > columnIndex) {
      assertion(*iter > 0);
      *iter -= 1;
      if (*iter == 0) {
        _matrixCols.erase(iter);
      }
      break;
    }
    iter++;
  }
}

void BaseQNPostProcessing::exportState(
    io::TXTWriter &writer)
{
}

void BaseQNPostProcessing::importState(
    io::TXTReader &reader)
{
}

int BaseQNPostProcessing::getDeletedColumns()
{
  return _nbDelCols;
}

int BaseQNPostProcessing::getLSSystemCols()
{
  int cols = 0;
  for (int col : _matrixCols) {
    cols += col;
  }
  if (_hasNodesOnInterface) {
    assertion(cols == _matrixV.cols(), cols, _matrixV.cols(), _matrixCols, _qrV.cols());
    assertion(cols == _matrixW.cols(), cols, _matrixW.cols());
  }

  return cols;
}

int BaseQNPostProcessing::getLSSystemRows()
{
  if (utils::MasterSlave::_masterMode || utils::MasterSlave::_slaveMode) {
    return _dimOffsets.back();
  }
  return _residuals.size();
}

void BaseQNPostProcessing::writeInfo(
    std::string s, bool allProcs)
{
  if (not utils::MasterSlave::_masterMode && not utils::MasterSlave::_slaveMode) {
    // serial post processing mode, server mode
    _infostringstream << s;

    // parallel post processing, master-slave mode
  } else {
    if (not allProcs) {
      if (utils::MasterSlave::_masterMode)
        _infostringstream << s;
    } else {
      _infostringstream << s;
    }
  }
  _infostringstream << std::flush;
}
}
}
} // namespace precice, cplscheme, impl
