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
#include <string.h>
//#include "utils/NumericalCompare.hpp"

#include <time.h>

namespace precice {
namespace cplscheme {
namespace impl {

tarch::logging::Log BaseQNPostProcessing::
      _log("precice::cplscheme::impl::BaseQNPostProcessing");

      
/* ----------------------------------------------------------------------------
 *     Constructor
 * ----------------------------------------------------------------------------
 */ 
BaseQNPostProcessing:: BaseQNPostProcessing
(
  double initialRelaxation,
  int    maxIterationsUsed,
  int    timestepsReused,
  double singularityLimit,
  std::vector<int> dataIDs,
  std::map<int,double> scalings)
:
  PostProcessing(),
  _initialRelaxation(initialRelaxation),
  _maxIterationsUsed(maxIterationsUsed),
  _timestepsReused(timestepsReused),
  _singularityLimit(singularityLimit),
  _dataIDs(dataIDs),
  _secondaryDataIDs(),
  _scalings(scalings),
  _firstIteration(true),
  _firstTimeStep(true),
  _hasNodesOnInterface(true),
  _oldXTilde(),
  //_secondaryOldXTildes(),
  _residuals(),
  _secondaryResiduals(),
  _scaledValues(),
  _scaledOldValues(),
  _oldResiduals(),
  _matrixV(),
  _matrixW(),
  _qrV(),
  _matrixVBackup(),
  _matrixWBackup(),
  _matrixColsBackup(),
  //_secondaryMatricesW(),
  _matrixCols(),
  _dimOffsets(),
  _infostream(),
  its(0),
  tSteps(0),
  deletedColumns(0)
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
   preciceCheck(tarch::la::greater(_singularityLimit, 0.0),
                "BaseQNPostProcessing()", "Singularity limit for QN "
                << "post-processing has to be larger than numerical zero ("
                << tarch::la::NUMERICAL_ZERO_DIFFERENCE << ")!");
   
   
  _infostream.open ("postProcessingInfo.txt", std::ios_base::out);
  _infostream << std::setprecision(16);
  _qrV.setfstream(&_infostream);
}


/* ----------------------------------------------------------------------------
 *     initialize
 * ----------------------------------------------------------------------------
 */ 
void BaseQNPostProcessing::initialize(DataMap& cplData) {
	preciceTrace1("initialize()", cplData.size());
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
	_scaledValues.append(DataValues(entries, init));
	_scaledOldValues.append(DataValues(entries, init));


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
		for (auto & elem : _dataIDs) {
			std::cout<<" Offsets:(vertex) \n"<<cplData[elem]->mesh->getVertexOffsets()<<std::endl;
		}
		for (size_t i = 0; i < _dimOffsets.size()-1; i++){
			int accumulatedNumberOfUnknowns = 0;
			for (auto & elem : _dataIDs) {
				auto & offsets = cplData[elem]->mesh->getVertexOffsets();
				accumulatedNumberOfUnknowns += offsets[i] * cplData[elem]->dimension;
			}
			_dimOffsets[i+1] = accumulatedNumberOfUnknowns;
		}

		// test that the computed number of unknown per proc equals the number of entries actually present on that proc
		int unknowns = _dimOffsets[utils::MasterSlave::_rank + 1] - _dimOffsets[utils::MasterSlave::_rank];
		assertion2(entries == unknowns, entries, unknowns);

		if(utils::MasterSlave::_masterMode){
			//ss<<" Offsets: \n"<<_dimOffsets<<std::endl;
			std::cout<<" Offsets:(unknowns) \n"<<_dimOffsets<<std::endl;
		}
		writeInfo(ss.str());
		ss.clear();

/*
		_dimOffsets.resize(utils::MasterSlave::_size + 1);
		if (utils::MasterSlave::_slaveMode) {
			utils::MasterSlave::_communication->send(((int) entries), 0);
			utils::MasterSlave::_communication->receive(&_dimOffsets[0], _dimOffsets.size(), 0);
		}
		if (utils::MasterSlave::_masterMode) {
			_dimOffsets[0] = 0;
			_dimOffsets[1] = entries;
			for (int rankSlave = 1; rankSlave < utils::MasterSlave::_size; rankSlave++) {
				int localDim = 0;
				utils::MasterSlave::_communication->receive(localDim, rankSlave);
				_dimOffsets[rankSlave + 1] = _dimOffsets[rankSlave] + localDim;
			}
			for (int rankSlave = 1; rankSlave < utils::MasterSlave::_size; rankSlave++) {
				utils::MasterSlave::_communication->send(&_dimOffsets[0], _dimOffsets.size(), rankSlave);
			}
			ss<<" Offsets (correct): \n"<<_dimOffsets<<std::endl;
			std::cout<<" Offsets (correct): \n"<<_dimOffsets<<std::endl;
		}
*/
	}

	// set the number of global rows in the QRFactorization. This is essential for the correctness in master-slave mode!
	_qrV.setGlobalRows(getLSSystemRows());


	// ---------------------------------------------------
	//debug output for master-slave mode
	if (utils::MasterSlave::_masterMode || utils::MasterSlave::_slaveMode) {
		ss << "processor [" << utils::MasterSlave::_rank<< "]: unknowns at interface: " << entries << std::endl;
		std::cout<<ss.str();
	}else{
		ss<< "unknowns at interface: " << entries << std::endl;
	}
	writeInfo(ss.str(), true);
	// ---------------------------------------------------

	// Fetch secondary data IDs, to be relaxed with same coefficients from IQN-ILS
	for (DataMap::value_type& pair : cplData){
		if (not utils::contained(pair.first, _dataIDs)) {
			_secondaryDataIDs.push_back(pair.first);
			int secondaryEntries = pair.second->values->size();
	//      _secondaryOldXTildes[pair.first].append(DataValues(secondaryEntries, init));
			_secondaryResiduals[pair.first].append(DataValues(secondaryEntries, init));
		}
	}

	// Append old value columns, if not done outside of post-processing already
	for (DataMap::value_type& pair : cplData){
		int cols = pair.second->oldValues.cols();
		if (cols < 1) { // Add only, if not already done
			//assertion1(pair.second->values->size() > 0, pair.first);
			pair.second->oldValues.append(
					CouplingData::DataMatrix(pair.second->values->size(), 1, 0.0));
		}
	}
}


/* ----------------------------------------------------------------------------
 *     scaling
 * ----------------------------------------------------------------------------
 */
void BaseQNPostProcessing:: scaling
(
  DataMap& cplData)
{
  preciceTrace("scaling()");

  int offset = 0;
  for (int id : _dataIDs){
    double factor = _scalings[id];
    preciceDebug("Scaling Factor " << factor << " for id: " << id);
    int size = cplData[id]->values->size();
    DataValues& values = *cplData[id]->values;
    DataValues& oldValues = cplData[id]->oldValues.column(0);
    for (int i=0; i < size; i++){
      _scaledValues[i+offset] = values[i]/factor;
      _scaledOldValues[i+offset] = oldValues[i]/factor;
    }
    offset += size;
  } 
}

/* ----------------------------------------------------------------------------
 *     undoScaling
 * ----------------------------------------------------------------------------
 */
void BaseQNPostProcessing:: undoScaling
(
  DataMap& cplData)
{
  preciceTrace("undoScaling()");
  
  int offset = 0;
  for (int id : _dataIDs){
    double factor = _scalings[id];
    int size = cplData[id]->values->size();
    preciceDebug("Copying values back, size: " << size);
    utils::DynVector& valuesPart = *(cplData[id]->values);
    utils::DynVector& oldValuesPart = cplData[id]->oldValues.column(0);
    for(int i=0; i < size; i++){
      valuesPart[i] = _scaledValues[i+offset]*factor;
      oldValuesPart[i] = _scaledOldValues[i+offset]*factor;
    }
    offset += size;
  }
}


/* ----------------------------------------------------------------------------
 *     updateDiffernceMatrices
 * ----------------------------------------------------------------------------
 */
void BaseQNPostProcessing::updateDifferenceMatrices(DataMap& cplData) {
	preciceTrace("updateDiffernceMatrices()");
	using namespace tarch::la;

	// Compute current residual: vertex-data - oldData
	_residuals = _scaledValues;
	_residuals -= _scaledOldValues;

	//if (_firstIteration && (_firstTimeStep || (_matrixCols.size() < 2))) {
	if (_firstIteration && _firstTimeStep){
		// do nothing: constant relaxation
	}else{
		preciceDebug("   Update Difference Matrices");
		if (not _firstIteration) {
			// Update matrices V, W with newest information

			assertion2(_matrixV.cols() == _matrixW.cols(), _matrixV.cols(), _matrixW.cols());
			assertion2(getLSSystemCols() <= _maxIterationsUsed,getLSSystemCols(), _maxIterationsUsed);
			
			if(2*getLSSystemCols() >= getLSSystemRows())
			  preciceWarning("updateDifferenceMatrices()", "The number of columns in the least squares system exceeded half the number of unknowns at the interface. The system will probably become bad or ill-conditioned and the quasi-Newton post processing may not converge. Maybe the number of allowed columns (maxIterationsUsed) should be limited.");

			DataValues deltaR = _residuals;
			deltaR -= _oldResiduals;

			DataValues deltaXTilde = _scaledValues;
			deltaXTilde -= _oldXTilde;

			bool columnLimitReached = getLSSystemCols() == _maxIterationsUsed;
			bool overdetermined = getLSSystemCols() <= getLSSystemRows();
			if (not columnLimitReached && overdetermined) {

				_matrixV.appendFront(deltaR);
				_matrixW.appendFront(deltaXTilde);

				// insert column deltaR = _residuals - _oldResiduals at pos. 0 (front) into the
				// QR decomposition and updae decomposition
				_qrV.pushFront(deltaR);

				_matrixCols.front()++;
			}else {
				_matrixV.shiftSetFirst(deltaR);
				_matrixW.shiftSetFirst(deltaXTilde);

				// inserts column deltaR at pos. 0 to the QR decomposition and deletes the last column
				// the QR decomposition of V is updated
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
		_oldXTilde = _scaledValues;   // Store x_tilde
	}

}



/* ----------------------------------------------------------------------------
 *     performPostProcessing
 * ----------------------------------------------------------------------------
 */
void BaseQNPostProcessing::performPostProcessing(DataMap& cplData) {
	preciceTrace2("performPostProcessing()", _dataIDs.size(), cplData.size());
	using namespace tarch::la;
	assertion2(_dataIDs.size() == _scalings.size(), _dataIDs.size(), _scalings.size());
	assertion2(_oldResiduals.size() == _oldXTilde.size(),_oldResiduals.size(), _oldXTilde.size());
	assertion2(_scaledValues.size() == _oldXTilde.size(),_scaledValues.size(), _oldXTilde.size());
	assertion2(_scaledOldValues.size() == _oldXTilde.size(),_scaledOldValues.size(), _oldXTilde.size());
	assertion2(_residuals.size() == _oldXTilde.size(),_residuals.size(), _oldXTilde.size());

	// scale data values (and secondary data values)
	scaling(cplData);

	/** update the difference matrices V,W  includes:
	 * scaling of values
	 * computation of residuals
	 * appending the difference matrices
	 */
	updateDifferenceMatrices(cplData);

	//if (_firstIteration && (_firstTimeStep || (_matrixCols.size() < 2))) {
	if (_firstIteration && _firstTimeStep){
		preciceDebug("   Performing underrelaxation");
		_oldXTilde = _scaledValues; // Store x tilde
		_oldResiduals = _residuals; // Store current residual

		// Perform constant relaxation
		// with residual: x_new = x_old + omega * res
		_residuals *= _initialRelaxation;
		_residuals += _scaledOldValues;
		_scaledValues = _residuals;

		computeUnderrelaxationSecondaryData(cplData);
	}else{
		preciceDebug("   Performing quasi-Newton Step");

		// If the previous time step converged within one single iteration, nothing was added
		// to the LS system matrices and they need to be restored from the backup at time T-2
		if (not _firstTimeStep && (getLSSystemCols() < 1) && (_timestepsReused == 0)) {
			preciceDebug("   Last time step converged after one iteration. Need to restore the matrices from backup.");

			_matrixCols = _matrixColsBackup;
			_matrixV = _matrixVBackup;
			_matrixW = _matrixWBackup;

			// re-computation of QR decomposition from _matrixV = _matrixVBackup
			// this occurs very rarely, to be precise, it occurs only if the coupling terminates
			// after the first iteration and the matrix data from time step t-2 has to be used
			_qrV.reset(_matrixV);
			// set the number of global rows in the QRFactorization. This is essential for the correctness in master-slave mode!
			_qrV.setGlobalRows(getLSSystemRows());
		}

		DataValues xUpdate(_residuals.size(), 0.0);

		/**
		 * compute quasi-Newton update
		 */
		computeQNUpdate(cplData, xUpdate);

		/**
		 * apply quasiNewton update
		 */
		_scaledValues = _scaledOldValues;  // = x^k
		_scaledValues += xUpdate;        // = x^k + delta_x
		_scaledValues += _residuals; // = x^k + delta_x + r^k

		// pending deletion: delete old V, W matrices if timestepsReused = 0
		// those were only needed for the first iteration (instead of underrelax.)
		if (_firstIteration && _timestepsReused == 0) {
			// save current matrix data in case the coupling for the next time step will terminate
			// after the first iteration (no new data, i.e., V = W = 0)
			if(getLSSystemCols() > 0){
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
	}

	// Undo scaling of data values and overwrite originals
	undoScaling(cplData);

	// number of iterations (usually equals number of columns in LS-system)
	its++;
	_firstIteration = false;
}

void BaseQNPostProcessing:: iterationsConverged
(
   DataMap & cplData)
{
  preciceTrace("iterationsConverged()");
  
  // debugging info, remove if not needed anymore:
  // -----------------------
  _infostream<<"\n ---------------- deletedColumns:"<<deletedColumns
		     <<"\n\n ### time step:"<<tSteps+1<<" ###"<<std::endl;
  its = 0;
  tSteps++;
  deletedColumns = 0;
  // -----------------------


  // writig l2 norm of converged configuration to info stream
  // -----------
  if(_firstTimeStep)
  {
    _infostream<<"l2-Norm of converged configuration after first time step:"<<std::endl;
    double l2norm = 0., oldl2norm = 0.;
    for (int id : _dataIDs)
    {
      l2norm = 0.; 
      double factor = _scalings[id];
    
      int size = cplData[id]->values->size();
      DataValues& values = *cplData[id]->values;
      DataValues unscaled(size, 0.0);
      for (int i=0; i < size; i++){
    	  unscaled(i) = values[i]/factor;
      }
      l2norm = utils::MasterSlave::l2norm(unscaled);
      if (id == _dataIDs[0]) oldl2norm = sqrt(l2norm);
      _infostream<<"  * "<<factor<<"  l2-norm: "<<sqrt(l2norm)<<"  of id: "<<id<<"\n"<<std::flush;
    } 
    _infostream<<"  * l2-norm ratio: "<<(double)oldl2norm/sqrt(l2norm)<<"\n"<<std::flush;
  }
  // -----------
  
  
  // the most recent differences for the V, W matrices have not been added so far
  // this has to be done in iterations converged, as PP won't be called any more if 
  // convergence was achieved
  scaling(cplData);
  updateDifferenceMatrices(cplData);
  undoScaling(cplData);

  
# ifdef Debug
  std::ostringstream stream;
  stream << "Matrix column counters: ";
  foreach (int cols, _matrixCols){
    stream << cols << ", ";
  }
  preciceDebug(stream.str());
# endif // Debug

  _firstTimeStep = false;
  if (_matrixCols.front() == 0){ // Did only one iteration
    _matrixCols.pop_front(); 
  }
  
  // doing specialized stuff for the corresponding post processing scheme after 
  // convergence of iteration i.e.:
  // - analogously to the V,W matrices, remove columns from matrices for secondary data
  // - save the old jacobian matrix
  specializedIterationsConverged(cplData);

  if (_timestepsReused == 0){
    
    /**
     * pending deletion (after first iteration of next time step
     * Using the matrices from the old time step for the first iteration
     * is better than doing underrelaxation as first iteration of every time step
     */
  }
  else if ((int)_matrixCols.size() > _timestepsReused){
    int toRemove = _matrixCols.back();
    assertion1(toRemove > 0, toRemove);
    preciceDebug("Removing " << toRemove << " cols from least-squares system with "
                 << getLSSystemCols() << " cols");
    assertion2(_matrixV.cols() == _matrixW.cols(), _matrixV.cols(), _matrixW.cols());
    assertion2(getLSSystemCols() > toRemove, getLSSystemCols(), toRemove);

    // remove columns
	for (int i=0; i < toRemove; i++){
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

void BaseQNPostProcessing:: exportState(io::TXTWriter& writer)
{
//  tarch::la::Vector<1,int> colSize(_matrixCols.size());
//  writer.write(colSize);
//  writer.write(_matrixCols);
//  writer.write(_oldXTilde);
//  writer.write(_oldResiduals);
//  writer.write(_matrixV);
//  writer.write(_matrixW);
}

void BaseQNPostProcessing:: importState(io::TXTReader& reader)
{
//  tarch::la::Vector<1,int> colSize;
//  reader.read(colSize);
//  _matrixCols.resize(colSize[0]);
//  reader.read(_oldXTilde);
//  reader.read(_oldResiduals);
//  for (int i=0; i<colSize[0]; i++ ){
//    // Using _oldXTilde to have a vector with appropriate size.
//    // Values are overwritten afterwards in file read.
//    _matrixV.append(_oldXTilde);
//    _matrixW.append(_oldXTilde);
//  }
//  reader.read(_matrixV);
//  reader.read(_matrixW);
//  if (colSize[0] > 1){
//    _firstIteration = false;
//  }
}

int BaseQNPostProcessing::getDeletedColumns()
{
	return deletedColumns;
}

int BaseQNPostProcessing::getLSSystemCols()
{
//	if(_hasNodesOnInterface){
//		return _matrixV.cols();
//	}
	int cols = 0;
	for (int col : _matrixCols){
		cols += col;
	}
	if(_hasNodesOnInterface){
		assertion4(cols == _matrixV.cols(), cols, _matrixV.cols(), _matrixCols, _qrV.cols());
		assertion2(cols == _matrixW.cols(), cols, _matrixW.cols());
	}

	return cols;
}

int BaseQNPostProcessing::getLSSystemRows()
{
	if(utils::MasterSlave::_masterMode || utils::MasterSlave::_slaveMode){
		return _dimOffsets.back();
	}
	return _residuals.size();
	//return _matrixV.rows();
}

void BaseQNPostProcessing:: removeMatrixColumn
(
  int columnIndex)
{
  preciceTrace2("removeMatrixColumn()", columnIndex, _matrixV.cols());


  // debugging information, can be removed
  deletedColumns++;

  assertion(_matrixV.cols() > 1);
  _matrixV.remove(columnIndex);
  _matrixW.remove(columnIndex);

  // remove corresponding column from dynamic QR-decomposition of _matrixV
  // Note: here, we need to delete the column, as we push empty columns
  //       for procs with no vertices in master-slave mode.
  _qrV.deleteColumn(columnIndex);

  
  // Reduce column count
  std::deque<int>::iterator iter = _matrixCols.begin();
  int cols = 0;
  while(iter != _matrixCols.end()) {
    cols += *iter;
    if(cols > columnIndex) {
      assertion(*iter > 0);
      *iter -= 1;
      if(*iter == 0) {
        _matrixCols.erase(iter);
      }
      break;
    }
    iter++;
  }
}

void BaseQNPostProcessing::writeInfo
(std::string s, bool allProcs)
{
	if(not utils::MasterSlave::_masterMode && not utils::MasterSlave::_slaveMode){
		// serial post processing mode, server mode
		_infostream<<s;

		// parallel post processing, master-slave mode
	}else{
		if(not allProcs){
			if(utils::MasterSlave::_masterMode) _infostream<<s;
		}else{
			_infostream<<s;
		}
	}
	_infostream<<std::flush;
}

}}} // namespace precice, cplscheme, impl
