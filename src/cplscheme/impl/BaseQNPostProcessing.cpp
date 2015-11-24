// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "BaseQNPostProcessing.hpp"
#include "cplscheme/CouplingData.hpp"
#include "utils/Globals.hpp"
#include "tarch/la/MatrixVectorOperations.h"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "utils/Dimensions.hpp"
#include "tarch/la/Scalar.h"
#include "io/TXTWriter.hpp"
#include "io/TXTReader.hpp"
#include "QRFactorization.hpp"
#include "utils/MasterSlave.hpp"
#include "utils/EventTimings.hpp"
#include <string.h>
#include <iostream>
#include <sstream>
//#include "utils/NumericalCompare.hpp"

#include <time.h>

using precice::utils::Event;

namespace precice {
namespace cplscheme {
namespace impl {

tarch::logging::Log BaseQNPostProcessing::
_log("precice::cplscheme::impl::BaseQNPostProcessing");


/* ----------------------------------------------------------------------------
 *     Constructor
 * ----------------------------------------------------------------------------
 */
BaseQNPostProcessing::BaseQNPostProcessing
(
    double initialRelaxation,
    bool forceInitialRelaxation,
    int maxIterationsUsed,
    int timestepsReused,
    int filter,
    double singularityLimit,
    std::vector<int> dataIDs,
    PtrPreconditioner preconditioner)
:
  PostProcessing(),
  _preconditioner(preconditioner),
  _initialRelaxation(initialRelaxation),
  _maxIterationsUsed(maxIterationsUsed),
  _timestepsReused(timestepsReused),
  _dataIDs(dataIDs),
  _secondaryDataIDs(),
  _firstIteration(true),
  _firstTimeStep(true),
  _hasNodesOnInterface(true),
  _forceInitialRelaxation(forceInitialRelaxation),
  _oldXTilde(),
  _residuals(),
  _secondaryResiduals(),
  _matrixV(),
  _matrixW(),
  _qrV(filter),
  _matrixCols(),
  _dimOffsets(),
  _values(),
  _oldValues(),
  _oldResiduals(),
  _filter(filter),
  _singularityLimit(singularityLimit),
  _designSpecification(),
  _matrixVBackup(),
  _matrixWBackup(),
  _matrixColsBackup(),
  its(0),
  tSteps(0),
  deletedColumns(0),
  _infostream()
{
  preciceCheck((_initialRelaxation > 0.0) && (_initialRelaxation <= 1.0),
      "BaseQNPostProcessing()",
      "Initial relaxation factor for QN post-processing has to "
      << "be larger than zero and smaller or equal than one!");
  preciceCheck(_maxIterationsUsed > 0, "BaseQNPostProcessing()",
      "Maximal iterations used for QN post-processing has to "
      << "be larger than zero!");
  preciceCheck(_timestepsReused >= 0, "BaseQNPostProcessing()",
      "Number of old timesteps to be reused for QN "
      << "post-processing has to be >= 0!");

  //_infostream.open("postProcessingInfo.txt", std::ios_base::out);
  //_infostream << std::setprecision(16);
  _qrV.setfstream(&_infostream);
}



/** ---------------------------------------------------------------------------------------------
 *         initialize()
 *
 * @brief: Initializes all the needed variables and data
 *  ---------------------------------------------------------------------------------------------
 */
void BaseQNPostProcessing::initialize(
    DataMap& cplData)
{
  preciceTrace1("initialize()", cplData.size());
  Event e("BaseQNPostProcessing::initialize", true, true); // time measurement, barrier

  size_t entries = 0;

  for (auto & elem : _dataIDs) {
    preciceCheck(utils::contained(elem, cplData), "initialize()",
        "Data with ID " << elem << " is not contained in data " "given at initialization!");
    entries += cplData[elem]->values->size();
  }

  _matrixCols.push_front(0);
  _firstIteration = true;
  _firstTimeStep = true;

  double init = 0.0;
  assertion(_oldXTilde.size() == 0);assertion(_oldResiduals.size() == 0);
  _oldXTilde.append(DataValues(entries, init));
  _oldResiduals.append(DataValues(entries, init));
  _residuals.append(DataValues(entries, init));
  _values.append(DataValues(entries, init));
  _oldValues.append(DataValues(entries, init));

  // if design specifiaction not initialized yet
  if (not (_designSpecification.size() > 0)) {
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
    //	std::cout<<" Offsets:(vertex) \n"<<cplData[elem]->mesh->getVertexOffsets()<<std::endl;
    //}
    for (size_t i = 0; i < _dimOffsets.size() - 1; i++) {
      int accumulatedNumberOfUnknowns = 0;
      for (auto & elem : _dataIDs) {
        auto & offsets = cplData[elem]->mesh->getVertexOffsets();
        accumulatedNumberOfUnknowns += offsets[i] * cplData[elem]->dimension;
      }
      _dimOffsets[i + 1] = accumulatedNumberOfUnknowns;
    }

    // test that the computed number of unknown per proc equals the number of entries actually present on that proc
    size_t unknowns = _dimOffsets[utils::MasterSlave::_rank + 1] - _dimOffsets[utils::MasterSlave::_rank];
    assertion2(entries == unknowns, entries, unknowns);
    writeInfo(ss.str());
    ss.clear();

  }

  // set the number of global rows in the QRFactorization. This is essential for the correctness in master-slave mode!
  _qrV.setGlobalRows(getLSSystemRows());

  // Fetch secondary data IDs, to be relaxed with same coefficients from IQN-ILS
  for (DataMap::value_type& pair : cplData) {
    if (not utils::contained(pair.first, _dataIDs)) {
      _secondaryDataIDs.push_back(pair.first);
      int secondaryEntries = pair.second->values->size();
      //      _secondaryOldXTildes[pair.first].append(DataValues(secondaryEntries, init));
      _secondaryResiduals[pair.first].append(DataValues(secondaryEntries, init));
    }
  }

  // Append old value columns, if not done outside of post-processing already
  for (DataMap::value_type& pair : cplData) {
    int cols = pair.second->oldValues.cols();
    if (cols < 1) { // Add only, if not already done
      //assertion1(pair.second->values->size() > 0, pair.first);
      pair.second->oldValues.append(
          CouplingData::DataMatrix(pair.second->values->size(), 1, 0.0));
    }
  }

  _preconditioner->initialize(entries);
}


/** ---------------------------------------------------------------------------------------------
 *         setDesignSpecification()
 *
 * @brief: sets a design specification for the fine model optimization problem
 *         i.e., x_star = argmin_x || f(x) - q ||
 *  ---------------------------------------------------------------------------------------------
 */
void BaseQNPostProcessing::setDesignSpecification
(
    Eigen::VectorXd& q)
{
  preciceTrace("setDesignSpecification()");
  assertion2(q.size() == _residuals.size(), q.size(), _residuals.size());
  _designSpecification = q;
}

/** ---------------------------------------------------------------------------------------------
 *         getDesignSpecification()
 *
 * @brief: Returns the design specification corresponding to the given coupling data.
 *         This information is needed for convergence measurements in the coupling scheme.
 *  ---------------------------------------------------------------------------------------------
 */        // TODO: change to call by ref when Eigen is used.
std::map<int, utils::DynVector> BaseQNPostProcessing::getDesignSpecification
(
  DataMap& cplData)
{
  std::map<int, utils::DynVector> designSpecifications;
  int off = 0;
  for (int id : _dataIDs) {
      int size = cplData[id]->values->size();
      utils::DynVector q(size, 0.0);
      for (int i = 0; i < size; i++) {
        q(i) = _designSpecification(i+off);
      }
      off += size;
      std::map<int, utils::DynVector>::value_type pair = std::make_pair(id, q);
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
void BaseQNPostProcessing::updateDifferenceMatrices
(
    DataMap& cplData)
{
  preciceTrace("updateDiffernceMatrices()");
  Event e(__func__, true, true); // time measurement, barrier
  using namespace tarch::la;

  // Compute current residual: vertex-data - oldData
  _residuals = _values;
  _residuals -= _oldValues;

  //if (_firstIteration && (_firstTimeStep || (_matrixCols.size() < 2))) {
  if (_firstIteration && (_firstTimeStep || _forceInitialRelaxation)) {
    // do nothing: constant relaxation
  } else {
    preciceDebug("   Update Difference Matrices");
    if (not _firstIteration) {
      // Update matrices V, W with newest information

      assertion2(_matrixV.cols() == _matrixW.cols(), _matrixV.cols(), _matrixW.cols());
      assertion2(getLSSystemCols() <= _maxIterationsUsed,getLSSystemCols(), _maxIterationsUsed);

      if (2 * getLSSystemCols() >= getLSSystemRows())
        preciceWarning("updateDifferenceMatrices()",
            "The number of columns in the least squares system exceeded half the number of unknowns at the interface. The system will probably become bad or ill-conditioned and the quasi-Newton post processing may not converge. Maybe the number of allowed columns (maxIterationsUsed) should be limited.");

      DataValues deltaR = _residuals;
      deltaR -= _oldResiduals;

      DataValues deltaXTilde = _values;
      deltaXTilde -= _oldXTilde;

      bool columnLimitReached = getLSSystemCols() == _maxIterationsUsed;
      bool overdetermined = getLSSystemCols() <= getLSSystemRows();
      if (not columnLimitReached && overdetermined) {

        _matrixV.appendFront(deltaR);
        _matrixW.appendFront(deltaXTilde);

        // insert column deltaR = _residuals - _oldResiduals at pos. 0 (front) into the
        // QR decomposition and updae decomposition

        //apply scaling here

        _preconditioner->apply(deltaR);
        _qrV.pushFront(deltaR);

        _matrixCols.front()++;
        }
      else {
        _matrixV.shiftSetFirst(deltaR);
        _matrixW.shiftSetFirst(deltaXTilde);

        // inserts column deltaR at pos. 0 to the QR decomposition and deletes the last column
        // the QR decomposition of V is updated
        _preconditioner->apply(deltaR);
        _qrV.pushFront(deltaR);
        _qrV.popBack();

        _matrixCols.front()++;
        _matrixCols.back()--;
        if(_matrixCols.back() == 0) {
          _matrixCols.pop_back();
        }
      }
    }
    _oldResiduals = _residuals;   // Store residuals
    _oldXTilde = _values;   // Store x_tilde
  }

}

/** ---------------------------------------------------------------------------------------------
 *         performPostProcessing()
 *
 * @brief: performs one iteration of the quasi Newton post processing. It steers the execution
 *         of fine and coarse model evaluations and also calls the coarse model optimization.
 *  ---------------------------------------------------------------------------------------------
 */
void BaseQNPostProcessing::performPostProcessing
(
    DataMap& cplData)
{
  preciceTrace2("performPostProcessing()", _dataIDs.size(), cplData.size());
  Event e(__func__, true, true); // time measurement, barrier

  using namespace tarch::la;
  assertion2(_oldResiduals.size() == _oldXTilde.size(),_oldResiduals.size(), _oldXTilde.size());
  assertion2(_values.size() == _oldXTilde.size(),_values.size(), _oldXTilde.size());
  assertion2(_oldValues.size() == _oldXTilde.size(),_oldValues.size(), _oldXTilde.size());
  assertion2(_residuals.size() == _oldXTilde.size(),_residuals.size(), _oldXTilde.size());

  // scale data values (and secondary data values)
  concatenateCouplingData(cplData);


  /** update the difference matrices V,W  includes:
   * scaling of values
   * computation of residuals
   * appending the difference matrices
   */
  updateDifferenceMatrices(cplData);

  if (_firstIteration && (_firstTimeStep || _forceInitialRelaxation)) {
    preciceDebug("   Performing underrelaxation");
    _oldXTilde = _values; // Store x tilde
    _oldResiduals = _residuals; // Store current residual

    // copying is removed when moving to Eigen
    DataValues q(_residuals.size(), 0.0);
    for (int i = 0; i < q.size(); i++)
      q(i) = _designSpecification(i) * _initialRelaxation;

    // Perform constant relaxation
    // with residual: x_new = x_old + omega * (res-q)
    _residuals *= _initialRelaxation;
    _residuals -= q;
    _residuals += _oldValues;
    _values = _residuals;

    computeUnderrelaxationSecondaryData(cplData);
  } else {
    preciceDebug("   Performing quasi-Newton Step");

    // If the previous time step converged within one single iteration, nothing was added
    // to the LS system matrices and they need to be restored from the backup at time T-2
    if (not _firstTimeStep && (getLSSystemCols() < 1) && (_timestepsReused == 0) && not _forceInitialRelaxation) {
      preciceDebug("   Last time step converged after one iteration. Need to restore the matrices from backup.");

      _matrixCols = _matrixColsBackup;
      _matrixV = _matrixVBackup;
      _matrixW = _matrixWBackup;

      // re-computation of QR decomposition from _matrixV = _matrixVBackup
      // this occurs very rarely, to be precise, it occurs only if the coupling terminates
      // after the first iteration and the matrix data from time step t-2 has to be used
      _qrV.reset(_matrixV, getLSSystemRows());
    }

    // subtract design specification from residuals, i.e., we want to minimize argmin_x|| r(x) - q ||
    assertion2(_residuals.size() == _designSpecification.size(), _residuals.size(), _designSpecification.size());
    for (int i = 0; i < _designSpecification.size(); i++)
          _residuals(i) -= _designSpecification(i);


    /**
     *  === update and apply preconditioner ===
     *
     * IQN-ILS would also work without W and xUpdate scaling, IQN-IMVJ unfortunately not
     * Note: here, the _residuals are H(x)- x - q, i.e., residual of the fixed-point iteration
     *       minus the design specification of the optimization problem (!= null if MM is used)
     */
    Event e_applyPrecond("applyPreconditioner", true, true); // time measurement, barrier
    _preconditioner->update(false, _values, _residuals);
    // TODO: evaluate whether the pure residual should be used for updating the preconditioner or residual - design specification
    _preconditioner->apply(_residuals);
    _preconditioner->apply(_matrixV);
    _preconditioner->apply(_matrixW);

    if(_preconditioner->requireNewQR()){
      if(not (_filter==PostProcessing::QR2FILTER)){ //for QR2 filter, there is no need to do this twice
        _qrV.reset(_matrixV, getLSSystemRows());
      }
      _preconditioner->newQRfulfilled();
    }
    e_applyPrecond.stop();

    // apply the configured filter to the LS system
    applyFilter();

    /**
     * compute quasi-Newton update
     */
    DataValues xUpdate(_residuals.size(), 0.0);
    computeQNUpdate(cplData, xUpdate);

    Event e_revertPrecond("revertPreconditioner", true, true); // time measurement, barrier
    _preconditioner->revert(xUpdate); //to compensate the W scaling
    _preconditioner->revert(_matrixW);
    _preconditioner->revert(_matrixV);
    _preconditioner->revert(_residuals);
    e_revertPrecond.stop();

    /**
     * apply quasiNewton update
     */
    _values = _oldValues;  // = x^k
    _values += xUpdate;        // = x^k + delta_x
    _values += _residuals; // = x^k + delta_x + r^k         note: residuals are _residuals - _designSpecifiaction at this point.
//    _values -= q; // = x^k + delta_x + r^k - q^k

    // TODO: maybe add design specification. Though, residuals are overwritten in the next iteration this would be a clearer and nicer code

    // pending deletion: delete old V, W matrices if timestepsReused = 0
    // those were only needed for the first iteration (instead of underrelax.)
    if (_firstIteration && _timestepsReused == 0 && not _forceInitialRelaxation) {
      // save current matrix data in case the coupling for the next time step will terminate
      // after the first iteration (no new data, i.e., V = W = 0)
      if (getLSSystemCols() > 0) {
        _matrixColsBackup = _matrixCols;
        _matrixVBackup = _matrixV;
        _matrixWBackup = _matrixW;
      }
      // if no time steps reused, the matrix data needs to be cleared as it was only needed for the
      // QN-step in the first iteration (idea: rather perform QN-step with information from last converged
      // time step instead of doing a underrelaxation)
      if (not _firstTimeStep) {
        _matrixV.clear();
        _matrixW.clear();
        _matrixCols.clear();
        _matrixCols.push_front(0); // vital after clear()
        _qrV.reset();
        // set the number of global rows in the QRFactorization. This is essential for the correctness in master-slave mode!
        _qrV.setGlobalRows(getLSSystemRows());
      }
    }

    if(std::isnan(utils::MasterSlave::l2norm(xUpdate))){
      preciceError(__func__, "The coupling iteration in time step "<<tSteps<<
          " failed to converge and NaN values occured throughout the coupling process. ");
    }
  }

  splitCouplingData(cplData);

  // number of iterations (usually equals number of columns in LS-system)
  its++;
  _firstIteration = false;
}


void BaseQNPostProcessing::applyFilter()
{
  preciceTrace1(__func__,_filter);
  Event e(__func__, true, true); // time measurement, barrier
  if (_filter == PostProcessing::NOFILTER) {
    // do nothing
  } else {
    // do: filtering of least-squares system to maintain good conditioning
    std::vector<int> delIndices(0);
    _qrV.applyFilter(_singularityLimit, delIndices, _matrixV);
    // start with largest index (as V,W matrices are shrinked and shifted
    for (int i = delIndices.size() - 1; i >= 0; i--) {

      removeMatrixColumn(delIndices[i]);

      std::stringstream ss;
      ss << "(updatedQR) removing linear dependent column " << delIndices[i] << "  time step: " << tSteps
          << " iteration: " << its << "\n" << std::endl;
      preciceDebug(ss.str());  writeInfo(ss.str());
    }
    assertion2(_matrixV.cols() == _qrV.cols(), _matrixV.cols(), _qrV.cols());
  }
}



void BaseQNPostProcessing::concatenateCouplingData
(
    DataMap& cplData)
{
  preciceTrace("concatenateCouplingData()");

  int offset = 0;
  for (int id : _dataIDs) {
    int size = cplData[id]->values->size();
    DataValues& values = *cplData[id]->values;
    DataValues& oldValues = cplData[id]->oldValues.column(0);
    for (int i = 0; i < size; i++) {
      _values[i + offset] = values[i];
      _oldValues[i + offset] = oldValues[i];
    }
    offset += size;
  }
}


void BaseQNPostProcessing::splitCouplingData
(
    DataMap& cplData)
{
  preciceTrace("splitCouplingData()");

  int offset = 0;
  for (int id : _dataIDs) {
    int size = cplData[id]->values->size();
    utils::DynVector& valuesPart = *(cplData[id]->values);
    utils::DynVector& oldValuesPart = cplData[id]->oldValues.column(0);
    for (int i = 0; i < size; i++) {
      valuesPart[i] = _values[i + offset];
      oldValuesPart[i] = _oldValues[i + offset];
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
void BaseQNPostProcessing::iterationsConverged
(
    DataMap & cplData)
{
  preciceTrace("iterationsConverged()");
  Event e(__func__, true, true); // time measurement, barrier

  // debugging info, remove if not needed anymore:
  // -----------------------
//  _infostream << "\n ---------------- deletedColumns:" << deletedColumns
//      << "\n\n ### time step:" << tSteps + 1 << " ###" << std::endl;
  its = 0;
  tSteps++;
  deletedColumns = 0;
  // -----------------------
 
  // the most recent differences for the V, W matrices have not been added so far
  // this has to be done in iterations converged, as PP won't be called any more if 
  // convergence was achieved
  concatenateCouplingData(cplData);
  updateDifferenceMatrices(cplData);

  Event e_upPrecond("iterConv: update Precond", true, true);
  // subtract design specification from residuals, i.e., we want to minimize argmin_x|| r(x) - q ||
  assertion2(_residuals.size() == _designSpecification.size(), _residuals.size(), _designSpecification.size());
  for (int i = 0; i < _designSpecification.size(); i++)
        _residuals(i) -= _designSpecification(i);

  _preconditioner->update(true, _values, _residuals);
  e_upPrecond.stop();

  // TODO: maybe add design specification. Though, residuals are overwritten in the next iteration this would be a clearer and nicer code
  _firstTimeStep = false;
  if (_matrixCols.front() == 0) { // Did only one iteration
    _matrixCols.pop_front();
  }

# ifdef Debug
  std::ostringstream stream;
  stream << "Matrix column counters: ";
  for (int cols: _matrixCols) {
    stream << cols << ", ";
  }
  preciceDebug(stream.str());
# endif // Debug

  // doing specialized stuff for the corresponding post processing scheme after 
  // convergence of iteration i.e.:
  // - analogously to the V,W matrices, remove columns from matrices for secondary data
  // - save the old jacobian matrix
  Event e_spec("specializedIterConverged", true, true);
  specializedIterationsConverged(cplData);
  e_spec.stop();
  
  if (_timestepsReused == 0) {

    if(_forceInitialRelaxation)
    {
      _matrixV.clear();
      _matrixW.clear();
      _qrV.reset();
      // set the number of global rows in the QRFactorization. This is essential for the correctness in master-slave mode!
      _qrV.setGlobalRows(getLSSystemRows());
      _matrixCols.clear(); // _matrixCols.push_front() at the end of the method.
    }else{
      /**
       * pending deletion (after first iteration of next time step
       * Using the matrices from the old time step for the first iteration
       * is better than doing underrelaxation as first iteration of every time step
       */
    }
  }
  else if ((int) _matrixCols.size() > _timestepsReused) {
    int toRemove = _matrixCols.back();
    assertion1(toRemove > 0, toRemove);
    preciceDebug("Removing " << toRemove << " cols from least-squares system with "
        << getLSSystemCols() << " cols");
    assertion2(_matrixV.cols() == _matrixW.cols(), _matrixV.cols(), _matrixW.cols());
    assertion2(getLSSystemCols() > toRemove, getLSSystemCols(), toRemove);

    // remove columns
    for (int i = 0; i < toRemove; i++) {
      _matrixV.remove(_matrixV.cols() - 1);
      _matrixW.remove(_matrixW.cols() - 1);
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
void BaseQNPostProcessing::removeMatrixColumn
(
    int columnIndex)
{
  preciceTrace2("removeMatrixColumn()", columnIndex, _matrixV.cols());

  // debugging information, can be removed
  deletedColumns++;

  assertion(_matrixV.cols() > 1);
  _matrixV.remove(columnIndex);
  _matrixW.remove(columnIndex);

  // Reduce column count
  std::deque<int>::iterator iter = _matrixCols.begin();
  int cols = 0;
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
    io::TXTWriter& writer)
{
}

void BaseQNPostProcessing::importState(
    io::TXTReader& reader)
{
}

int BaseQNPostProcessing::getDeletedColumns()
{
  return deletedColumns;
}

int BaseQNPostProcessing::getLSSystemCols()
{
  int cols = 0;
  for (int col : _matrixCols) {
    cols += col;
  }
  if (_hasNodesOnInterface) {
    assertion4(cols == _matrixV.cols(), cols, _matrixV.cols(), _matrixCols, _qrV.cols());
    assertion2(cols == _matrixW.cols(), cols, _matrixW.cols());
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


void BaseQNPostProcessing::writeInfo
(
    std::string s, bool allProcs)
{
  if (not utils::MasterSlave::_masterMode && not utils::MasterSlave::_slaveMode) {
    // serial post processing mode, server mode
//    _infostream << s;

    // parallel post processing, master-slave mode
  } else {
    if (not allProcs) {
//      if (utils::MasterSlave::_masterMode) _infostream << s;
    } else {
//      _infostream << s;
    }
  }
//  _infostream << std::flush;
}

}}} // namespace precice, cplscheme, impl
