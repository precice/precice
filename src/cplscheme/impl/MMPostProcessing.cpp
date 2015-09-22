/*
 * MMPostProcessing.cpp
 *
 *  Created on: Sep 18, 2015
 *      Author: Klaudius Scheufele
 */

// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "MMPostProcessing.hpp"
#include "cplscheme/CouplingData.hpp"
#include "utils/Globals.hpp"
//#include "tarch/la/MatrixVectorOperations.h"
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

tarch::logging::Log MMPostProcessing::
      _log("precice::cplscheme::impl::MMPostProcessing");

/* ----------------------------------------------------------------------------
 *     Constructor
 * ----------------------------------------------------------------------------
 */
MMPostProcessing:: MMPostProcessing
(
   impl::PtrPostProcessing coarseModelOptimization,
   int    maxIterationsUsed,
   int    timestepsReused,
   int    filter,
   double singularityLimit,
   bool estimateJacobian,
   std::vector<int>    fineDataIDs,
   std::vector<int>    coarseDataIDs,
   std::map<int,double>    scalings)
:
  PostProcessing(),
  _coarseModelOptimization(coarseModelOptimization),
  _maxIterationsUsed(maxIterationsUsed),
  _timestepsReused(timestepsReused),
  _singularityLimit(singularityLimit),
  _designSpecification(),
  _coarseModel_designSpecification(),
  _nextModelToEvaluate(),
  _fineDataIDs(fineDataIDs),
  _coarseDataIDs(coarseDataIDs),
  _dataIDs(),
  _secondaryDataIDs(),
  _scalings(scalings),
  _firstIteration(true),
  _firstTimeStep(true),
  _estimateJacobian(estimateJacobian),
  _fineResiduals(),
  _coarseResiduals(),
  _fineOldResiduals(),
  _coarseOldResiduals(),
  _outputFineModelScaled(),
  _outputCoarseModelScaled(),
  _input_Xstar(),
  _matrixF(),
  _matrixC(),
  _MMMappingMatrix(),
  _MMMappingMatrix_prev(),
  _matrixCols(),
  _dimOffsets(),
  its(0),
  tSteps(0),
  deletedColumns(0),
  _filter(filter)
{
   preciceCheck(_maxIterationsUsed > 0, "BaseQNPostProcessing()",
                "Maximal iterations used for QN post-processing has to "
                << "be larger than zero!");
   preciceCheck(_timestepsReused >= 0, "BaseQNPostProcessing()",
                "Number of old timesteps to be reused for QN "
                 << "post-processing has to be >= 0!");
}


/* ----------------------------------------------------------------------------
 *     initialize
 * ----------------------------------------------------------------------------
 */
void MMPostProcessing::initialize(DataMap& cplData) {
	preciceTrace1("initialize()", cplData.size());
	size_t entries = 0;
	size_t coarseEntries = 0;

	assertion2(_fineDataIDs.size() == _coarseDataIDs.size(), _fineDataIDs.size(), _coarseDataIDs.size());
	assertion1(_dataIDs.size() == 0, _dataIDs.size());

	_dataIDs.insert(_dataIDs.end(), _fineDataIDs.begin(), _fineDataIDs.end());
	_dataIDs.insert(_dataIDs.end(), _coarseDataIDs.begin(), _coarseDataIDs.end());


	for (auto & elem : _fineDataIDs) {
		preciceCheck(utils::contained(elem, cplData), "initialize()",
				"Data with ID " << elem << " is not contained in data " "given at initialization!");
		entries += cplData[elem]->values->size();
	}

	for (auto & elem : _coarseDataIDs) {
		preciceCheck(utils::contained(elem, cplData), "initialize()",
				"Data with ID " << elem << " is not contained in data " "given at initialization!");
		coarseEntries += cplData[elem]->values->size();
	}

	// the coarse model also uses the fine mesh (only evaluation in solver is on coarse model)
	assertion2(entries == coarseEntries, entries, coarseEntries);

	_matrixCols.push_front(0);
	_firstIteration = true;
	_firstTimeStep = true;

	double init = 0.0;
	assertion(_coarseOldResiduals.size() == 0);assertion(_fineOldResiduals.size() == 0);
	_coarseOldResiduals = Eigen::VectorXd::Zero(entries);
	_fineOldResiduals = Eigen::VectorXd::Zero(entries);
	_fineResiduals = Eigen::VectorXd::Zero(entries);
	_coarseResiduals = Eigen::VectorXd::Zero(entries);
	_outputFineModelScaled = Eigen::VectorXd::Zero(entries);
	_outputCoarseModelScaled = Eigen::VectorXd::Zero(entries);
	_input_Xstar = Eigen::VectorXd::Zero(entries);

	// if design specifiaction not initialized yet
	if(not (_designSpecification.size() > 0)){
		_designSpecification = Eigen::VectorXd::Zero(entries);
	}
	_coarseModel_designSpecification = Eigen::VectorXd::Zero(entries);

	/**
	 *  make dimensions public to all procs,
	 *  last entry _dimOffsets[MasterSlave::_size] holds the global dimension, global,n
	 */
	if (utils::MasterSlave::_masterMode || utils::MasterSlave::_slaveMode) {
		assertion(utils::MasterSlave::_communication.get() != NULL);
		assertion(utils::MasterSlave::_communication->isConnected());

		/** provide vertex offset information for all processors
		 *  mesh->getVertexOffsets() provides an array that stores the number of mesh vertices on each processor
		 *  This information needs to be gathered for all meshes. To get the number of respective unknowns of a specific processor
		 *  we need to multiply the number of vertices with the dimensionality of the vector-valued data for each coupling data.
		 */
		_dimOffsets.resize(utils::MasterSlave::_size + 1);
		_dimOffsets[0] = 0;
		for (size_t i = 0; i < _dimOffsets.size()-1; i++){
			int accumulatedNumberOfUnknowns = 0;
			for (auto & elem : _fineDataIDs) {
				auto & offsets = cplData[elem]->mesh->getVertexOffsets();
				accumulatedNumberOfUnknowns += offsets[i] * cplData[elem]->dimension;
			}
			_dimOffsets[i+1] = accumulatedNumberOfUnknowns;
		}

		// test that the computed number of unknown per proc equals the number of entries actually present on that proc
		size_t unknowns = _dimOffsets[utils::MasterSlave::_rank + 1] - _dimOffsets[utils::MasterSlave::_rank];
		assertion2(entries == unknowns, entries, unknowns);
	}

	if(_estimateJacobian){
		_MMMappingMatrix = Eigen::MatrixXd::Zero(getLSSystemRows(), entries);
		// do not initialize Tkprev (MMMappingMatrix_prev), as we need a different constructing rule
		// in the first step (see computeCoarseModelDesignSpecifiaction).

		//_MMMappingMatrix_prev = Eigen::MatrixXd::Identity(getLSSystemRows(), entries);
	}


	/**
	// Fetch secondary data IDs, to be relaxed with same coefficients from IQN-ILS
	for (DataMap::value_type& pair : cplData){
		if (not utils::contained(pair.first, _dataIDs)) {
			_secondaryDataIDs.push_back(pair.first);
			int secondaryEntries = pair.second->values->size();
	//      _secondaryOldXTildes[pair.first] Eigen::VectorXd::Zero(secondaryEntries);
			_secondaryResiduals[pair.first] = Eigen::VectorXd::Zero(secondaryEntries);
		}
	}
	*/

	// Append old value columns, if not done outside of post-processing already
	for (DataMap::value_type& pair : cplData){
		int cols = pair.second->oldValues.cols();
		if (cols < 1) { // Add only, if not already done
			//assertion1(pair.second->values->size() > 0, pair.first);
			pair.second->oldValues.append(CouplingData::DataMatrix(pair.second->values->size(), 1, 0.0));
		}
	}
}


/* ----------------------------------------------------------------------------
 *     scaling
 * ----------------------------------------------------------------------------
 */
void MMPostProcessing:: scaling
(
  DataMap& cplData)
{
  preciceTrace("scaling()");

  // scale design specification according to scaling of corresponding data
  bool isSet_designSpec = ((_designSpecification.size() > 0) && (_designSpecification.norm() > 1.0e-15 ));
  if(isSet_designSpec) assertion2(_outputFineModelScaled.size() == _designSpecification.size(), _outputFineModelScaled.size(), _designSpecification.size());


  int offset = 0;
  for (int id : _fineDataIDs){
    double factor = _scalings[id];
    preciceDebug("Scaling Factor " << factor << " for id: " << id);
    int size = cplData[id]->values->size();
    utils::DynVector& values = *cplData[id]->values;
    utils::DynVector& oldValues = cplData[id]->oldValues.column(0);
    for (int i=0; i < size; i++){
      _outputFineModelScaled[i+offset] = values[i]/factor;
      if(isSet_designSpec) _designSpecification[i+offset] = _designSpecification[i+offset]/factor;
      // ignore input from fine model as it must be exactly the
      // same as the input for the coarse model, if the fine model is evaluated

      //_scaledOldValues[i+offset] = oldValues[i]/factor;
    }
    offset += size;
  }
  offset = 0;
  for (int id : _coarseDataIDs){
    double factor = _scalings[id];
    preciceDebug("Scaling Factor " << factor << " for id: " << id);
    int size = cplData[id]->values->size();
    utils::DynVector& values = *cplData[id]->values;
    utils::DynVector& oldValues = cplData[id]->oldValues.column(0);
    for (int i=0; i < size; i++){
      _outputCoarseModelScaled[i+offset] = values[i]/factor;
	  _input_Xstar[i+offset] = oldValues[i]/factor;
    }
    offset += size;
  }
}

/* ----------------------------------------------------------------------------
 *     undoScaling (not needed in MM)
 * ----------------------------------------------------------------------------
 */
/*
void MMPostProcessing:: undoScaling
(
  DataMap& cplData)
{
  preciceTrace("undoScaling()");

  int offset = 0;
  for (int id : _fineDataIDs){
    double factor = _scalings[id];
    int size = cplData[id]->values->size();
    preciceDebug("Copying values back, size: " << size);
    utils::DynVector& valuesPart = *(cplData[id]->values);
    utils::DynVector& oldValuesPart = cplData[id]->oldValues.column(0);
    for(int i=0; i < size; i++){
      // write new coarse model solution back as input data for the fine model evaluation
      // _input_xStar needs to be updated in each iteration
      valuesPart[i] = _input_Xstar[i+offset]*factor;
     // oldValuesPart[i] = _scaledOldValues[i+offset]*factor; // not needed, as not modified
    }
    offset += size;
  }
  offset = 0;
  for (int id : _coarseDataIDs){
    double factor = _scalings[id];
    int size = cplData[id]->values->size();
    preciceDebug("Copying values back, size: " << size);
    utils::DynVector& valuesPart = *(cplData[id]->values);
    utils::DynVector& oldValuesPart = cplData[id]->oldValues.column(0);
    for(int i=0; i < size; i++){
      // write new coarse model solution back as input data for the fine model evaluation
      // _input_xStar needs to be updated in each iteration
      valuesPart[i] = _input_Xstar[i+offset]*factor;
      // oldValuesPart[i] = _coarseScaledOldValues[i+offset]*factor; // not needed, as not modified
    }
    offset += size;
  }
}
*/

/* ----------------------------------------------------------------------------
 *     registerSolutionCoarseModelOptimization
 * ----------------------------------------------------------------------------
 */
void MMPostProcessing::registerSolutionCoarseModelOptimization
(
   DataMap& cplData)
{
	// extract new solution x_star from coarse model optimization problem from coarse cplData
	int off = 0;
	for (int id : _coarseDataIDs){
	  int size = cplData[id]->values->size();
	  utils::DynVector& valuesPart = *(cplData[id]->values);
	  for(int i=0; i < size; i++){
		   _input_Xstar[i+off] = valuesPart[i];
	   }
	   off += size;
	}
	// register new solution x_star from coarse model optimization problem as input
	// to fine model evaluation.
	off = 0;
	for (int id : _fineDataIDs){
	   int size = cplData[id]->values->size();
	   utils::DynVector& valuesPart = *(cplData[id]->values);
	   for(int i=0; i < size; i++){
	     // write new coarse model solution back as input data for the fine model evaluation
	     // _input_xStar needs to be updated in each iteration
	     valuesPart[i] = _input_Xstar[i+off];
	   }
	   off += size;
	}
}

void MMPostProcessing::setDesignSpecification(Eigen::VectorXd& q)
{
	preciceTrace("setDesignSpecification()");
	assertion2(q.size() == _residuals.size(), q.size(), _residuals.size());
	_designSpecification = (q.norm() <= 1.0e-15) ? Eigen::VectorXd::Zero(_residuals.size()) : q;

	// only in the first step, the coarse model design specification equals the design specification
	// for the overall objective function (initial coarse solution)
	if(_firstTimeStep) _coarseModel_designSpecification = _designSpecification;
}


/** ----------------------------------------------------------------------------
 *     updateDiffernceMatrices
 *  ----------------------------------------------------------------------------
 */
void MMPostProcessing::updateDifferenceMatrices(DataMap& cplData) {

   preciceTrace("updateDiffernceMatrices()");

   /**
    * Compute current residual: vertex-data - oldData
    */
   _fineResiduals = _scaledValues - _input_Xstar;
   _coarseResiduals -= _coarseScaledValues - _input_Xstar;

   /**
    * Update matrices C, F with newest information
    */
   if(not _firstIteration)
   {
      preciceDebug("   Update Difference Matrices C and F with coarse and fine model responses");
      assertion2(_matrixF.cols() == _matrixC.cols(), _matrixF.cols(), _matrixC.cols());
      assertion2(getLSSystemCols() <= _maxIterationsUsed,getLSSystemCols(), _maxIterationsUsed);

      if(2*getLSSystemCols() >= getLSSystemRows())
	     preciceWarning("updateDifferenceMatrices()", "The number of columns in the least squares system exceeded half the number of unknowns at the interface. The system will probably become bad or ill-conditioned and the quasi-Newton post processing may not converge. Maybe the number of allowed columns (maxIterationsUsed) should be limited.");

      Eigen::VectorXd colF = _fineResiduals - _fineOldResiduals;
      Eigen::VectorXd colC = _coarseResiduals - _coarseOldResiduals;

      bool columnLimitReached = getLSSystemCols() == _maxIterationsUsed;
      bool overdetermined = getLSSystemCols() <= getLSSystemRows();
      if (not columnLimitReached && overdetermined) {

	     appendFront(_matrixF, colF);
	     appendFront(_matrixC, colC);

	     _matrixCols.front()++;
      }else {
	     shiftSetFirst(_matrixF, colF);
	     shiftSetFirst(_matrixC, colC);

	     _matrixCols.front()++;
	     _matrixCols.back()--;
	     if(_matrixCols.back() == 0) {
		    _matrixCols.pop_back();
	     }
      }
   }

   /**
    *  Store residuals
    */
   _fineOldResiduals = _fineResiduals;
   _coarseOldResiduals = _coarseResiduals;

}



/* ----------------------------------------------------------------------------
 *     performPostProcessing
 * ----------------------------------------------------------------------------
 */
void MMPostProcessing::performPostProcessing(DataMap& cplData) {
	preciceTrace2("performPostProcessing()", _dataIDs.size(), cplData.size());

	assertion2(_dataIDs.size() == _scalings.size(), _dataIDs.size(), _scalings.size());
	assertion2(_fineOldResiduals.size() == _fineResiduals.size(),_oldResiduals.size(), _fineResiduals.size());
	assertion2(_coarseResiduals.size() == _fineResiduals.size(),_scaledValues.size(), _fineResiduals.size());
	assertion2(_coarseOldResiduals.size() == _fineResiduals.size(),_scaledOldValues.size(), _fineResiduals.size());
	assertion2(_scaledValues.size() == _fineResiduals.size(),_residuals.size(), _fineResiduals.size());
	assertion2(_scaledOldValues.size() == _fineResiduals.size(),_residuals.size(), _fineResiduals.size());


	if(_nextModelToEvaluate == BaseCouplingScheme::ModelResolution::fineModel){

		/**
		 * assume the coarse model and the fine model has been evaluated for the new coarse model
		 * solution _input_Xstar, obtained in the coarse model optimization step.
		 */

		// scale data values (and secondary data values)
		scaling(cplData);

		// update the difference matrices with the newest residual deltas
		updateDifferenceMatrices(cplData);

		/** compute the new design specification for the coarse model optimization
		 *  updates: _coarseModel_designSpecification
		 *           i.e., qk = c(xk) - Tk * ( f(xk) - q )
		 *  updates: _MMMappingMatrix, i.e.,
		 *           Tk = Tkprev + (C - Tkprev * F) * pseudoInv_F,        iff Tkprev exists
		 *           Tk = C * pseudoInv_F + (I - Uc*Uc^T)*(I - Uf*Uf^T),  else (in the first step or after rescaling)
		 */
		computeCoarseModelDesignSpecifiaction();

		/**
		 * now, the difference matrices for the MM mapping as well as the design specification for the coarse
		 * model optimization problem are updated (also Jacobian of MM mapping matrix if required).
		 * next step: coarse model optimization, set the steering variable accordingly
		 */
		_nextModelToEvaluate = BaseCouplingScheme::ModelResolution::coarseModel;

		// Undo scaling of data values and overwrite originals
		undoScaling(cplData);

		// one MM iteration completed
	    its++;
		_firstIteration = false;
	}else{
		// view on coarse coupling data only
		DataMap coarseCplData;
		for (int id : _coarseDataIDs){
		  DataMap::value_type pair = std::make_pair (id, cplData[id]);
		  coarseCplData.insert(pair);
		}

		// NO SCALING of input/output data. This is done in the coarse optimization routine, only for the coarse cplData.
		// Though, the coarse model design specification is computed with scaled data and needs to be re-scaled to normal.
		// It is to be scaled again in the coarse model optimization scheme.
		if(not (_firstIteration && _firstTimeStep)) // not in very first iteration where we solve for a initial coarse solution with unscaled design specification
		{
		  int offset = 0;
	      for (int id : _fineDataIDs){
		    double factor = _scalings[id];
		    for (int i=0; i < cplData[id]->values->size(); i++){
			    _coarseModel_designSpecification[i+offset] = _coarseModel_designSpecification[i+offset]*factor;
		    }
		    offset += cplData[id]->values->size();
	      }
		}

		/** perform the coarse model optimization, determine x_star
		 *          x_k+1 = argmin_x || c(x) - q_k ||
		 *        ------------------------------------
		 *  In the first time step, the coarse model design specification is equal to the
		 *  design specification of the overall objective function. Here, a initial coarse
		 *  model solution is obtained for a initial guess.
		 */
		_coarseModelOptimization->optimize(coarseCplData, _coarseModel_designSpecification);

		/**
		 * after the coarse model optimization has converged successfully,
		 * the fine and the coarse model has to be evaluated for the new solution x_star.
		 * Hence, x_star needs to be copied to the fine model input values.
		 */
		registerSolutionCoarseModelOptimization(cplData);
	}
}


/** @brief: computes/updates the design specification for the coarse model optimization problem
 *     	   i. e., q_k = c(x_k) - T_k * (f(x_k) - q), q = 0 is the fine model design specification
 */
void MMPostProcessing::computeCoarseModelDesignSpecifiaction()
{

	/** update design specification
	 *  alpha = (f(x) - q),
	 *  q_k   = c(x)
	 */
	Eigen::VectorXd alpha = _fineResiduals - _designSpecification;
	_coarseModel_designSpecification = _coarseResiduals;


	// if residual differences are available for fine and coarse model
	// (either from previous iterations or from previous time steps or both)
	if(getLSSystemCols() > 0)
	{
		// compute SVDs of _matrixF and _matriC
		Eigen::VectorXd S_F, S_C;
		Eigen::MatrixXd V_F, U_F, V_C, U_C, Sigma_F, pseudoSigma_F;

		// Remove dependent columns of _matrixC and _matrixF
		int nbRemoveCols = 1;
		while ( nbRemoveCols > 0 )
		{
			nbRemoveCols = 0;
			if ( getLSSystemCols() == 0 )
				break;

			// Calculate singular value decomposition with Eigen
			Eigen::JacobiSVD<matrix> svd( _matrixF, Eigen::ComputeThinU | Eigen::ComputeThinV );
			Eigen::VectorXd singularValues = svd.singularValues();

			for ( int i = 0; i < singularValues.rows(); i++ ){
				if ( std::abs( singularValues( i ) ) <= singularityLimit ){

					// Remove the column from _matrixC and _matrixF
					removeMatrixColumn(i - nbRemoveCols);
					nbRemoveCols++;
				}
			}
			if ( nbRemoveCols )
				preciceDebug("Manifold mapping: remove " << nbRemoveCols << " columns from the Jacobian matrices");
		}

		assert( _matrixF.cols() == _matrixC.cols() );

		if ( getLSSystemCols() > 0 )
		{
			// Calculate singular value decomposition with Eigen
			Eigen::JacobiSVD<matrix> svd_C( _matrixC, Eigen::ComputeThinU | Eigen::ComputeThinV );
			Eigen::JacobiSVD<matrix> svd_F( _matrixF, Eigen::ComputeThinU | Eigen::ComputeThinV );

			Eigen::MatrixXd pseudoSigma_F = svd_F.singularValues().asDiagonal();

			for ( int i = 0; i < pseudoSigma_F.cols(); i++ )
				pseudoSigma_F( i, i ) = 1.0 / pseudoSigma_F( i, i );

			Eigen::MatrixXd pseudoMatrixF = svd_F.matrixV() * pseudoSigma_F * svd_F.matrixU().transpose();

			U_F = svd_F.matrixU();
			U_C = svd_C.matrixU();

			// if Jacobian matrix of MM mapping matrix is not set up explicitly, perform
			// matrix-vector update
			if (not _estimateJacobian)
			{
				Eigen::VectorXd beta = U_F * (U_F.transpose() * alpha);

				_coarseModel_designSpecification -= alpha;
				_coarseModel_designSpecification -= _matrixC * (pseudoMatrixF * alpha);
				_coarseModel_designSpecification += U_C * ( U_C.transpose() * (alpha - beta) );
				_coarseModel_designSpecification += beta;
			}

			// Jacobian matrix of MM mapping matrix is estimated and set up explicitly
			// multi-vector method for update of Jacobian with implicit incorporation of
			// information from previous time steps.
			if (_estimateJacobian)
			{
				Eigen::MatrixXd I = Eigen::MatrixXd::Identity( _matrixF.rows(), _matrixF.rows() );

				// if previous Jacobian exists, i.e., no re-scaling and not first estimation
				if ( _MMMappingMatrix_prev.rows() == getLSSystemRows() ){
					_MMMappingMatrix = _MMMappingMatrix_prev + (_matrixC - _MMMappingMatrix_prev * _matrixF) * pseudoMatrixF;

				// if no previous Jacobian exists, set up Jacobian with IQN-ILS update rule + stabilization term
				}else{
					_MMMappingMatrix = _matrixC * pseudoMatrixF + ( I - U_C * U_C.transpose() ) * ( I - U_F * U_F.transpose() );
			    }

				// compute new design specification for coarse model optimization: qk = c(x) - Tk( f(x) - q )
				_coarseModel_designSpecification -= _MMMappingMatrix * alpha;
			}
		}
	}

	// if no residual differences for the fine and coarse model are given so far
	if((_firstIteration && _firstTimeStep) || getLSSystemCols() <= 0)
	{
	   assertion1(getLSSystemCols() <= 0, getLSSystemCols());

	   if(_estimateJacobian && (_MMMappingMatrix_prev.rows() == getLSSystemRows()))
	   {
		   _coarseModel_designSpecification -= _MMMappingMatrix_prev * alpha;
	   }else{
		   _coarseModel_designSpecification -= alpha;
	   }
	}

}



void MMPostProcessing:: iterationsConverged
(
   DataMap & cplData)
{
  preciceTrace("iterationsConverged()");

  its = 0;
  tSteps++;
  deletedColumns = 0;

  // the most recent differences for the F, C matrices have not been added so far
  // this has to be done in iterations converged, as PP won't be called any more if
  // convergence was achieved
  scaling(cplData);
  updateDifferenceMatrices(cplData);
  // no undoScaling() needed, as input/output data is not modified


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

  if (_timestepsReused == 0){
    _matrixF.resize(0,0);
    _matrixC.resize(0,0);
    _matrixCols.clear();
  }
  else if ((int)_matrixCols.size() > _timestepsReused){
    int toRemove = _matrixCols.back();
    assertion1(toRemove > 0, toRemove);
    preciceDebug("Removing " << toRemove << " cols from mannifold mapping least-squares system with "<< getLSSystemCols() << " cols");
    assertion2(_matrixF.cols() == _matrixC.cols(), _matrixF.cols(), _matrixC.cols());
    assertion2(getLSSystemCols() > toRemove, getLSSystemCols(), toRemove);

    // remove columns
	for (int i=0; i < toRemove; i++){
		removeColumnFromMatrix(_matrixF, _matrixF.cols() - 1);
		removeColumnFromMatrix(_matrixC, _matrixC.cols() - 1);
    }
    _matrixCols.pop_back();
  }

  _matrixCols.push_front(0);
  _firstIteration = true;
}

void MMPostProcessing:: exportState(io::TXTWriter& writer)
{
}

void MMPostProcessing:: importState(io::TXTReader& reader)
{
}

int BaseQNPostProcessing::getDeletedColumns()
{
	return deletedColumns;
}

int MMPostProcessing::getLSSystemCols()
{
	int cols = 0;
	for (int col : _matrixCols){
		cols += col;
	}
	//if(_hasNodesOnInterface){
	//	assertion3(cols == _matrixF.cols(), cols, _matrixF.cols(), _matrixCols);
	//	assertion2(cols == _matrixC.cols(), cols, _matrixC.cols());
	//}

	return cols;
}

int MMPostProcessing::getLSSystemRows()
{
	if(utils::MasterSlave::_masterMode || utils::MasterSlave::_slaveMode){
		return _dimOffsets.back();
	}
	return _fineResiduals.size();
	//return _matrixF.rows();
}


void MMPostProcessing:: removeMatrixColumn
(
  int columnIndex)
{
  preciceTrace2("removeMatrixColumn()", columnIndex, _matrixF.cols());


  // debugging information, can be removed
  deletedColumns++;

  assertion(_matrixF.cols() > 1);
  removeColumnFromMatrix(_matrixF, columnIndex);
  removeColumnFromMatrix(_matrixC, columnIndex);

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


void MMPostProcessing::shiftSetFirst
(
   Eigen::MatrixXd& A, Eigen::VectorXd& v)
{
	assertion2(v.size() == A.rows(), v.size(), A.rows());
	int n = A.rows(), m = A.cols();
	A.bottomRightCorner(n,m-1) = A.topLeftCorner(n,m-1);
	A.col(0) = v;
}

void MMPostProcessing::appendFront
(
   Eigen::MatrixXd& A, Eigen::VectorXd& v)
{
   int n = A.rows(), m = A.cols();
   if(n <= 0 && m <= 0){
	   A = v;
   }else{
	  assertion2(v.size() == n, v.size(), A.rows());
      A.conservativeResize(n,m+1);
      A.bottomRightCorner(n,m) = A.topLeftCorner(n,m);
      A.col(0) = v;
   }
}

void MMPostProcessing::removeColumnFromMatrix
(
   Eigen::MatrixXd& A, int col)
{
	assertion2(col < A.cols() && col >= 0, col, A.cols())
    for ( int j = col; j < A.cols() - 1; j++ )
        A.col( j ) = A.col( j + 1 );

    A.conservativeResize( A.rows(), A.cols() - 1 );
}



}}} // namespace precice, cplscheme, impl



