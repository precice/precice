// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "IQNILSPostProcessing.hpp"
#include "cplscheme/CouplingData.hpp"
#include "utils/Globals.hpp"
#include "tarch/la/GramSchmidt.h"
#include "tarch/la/MatrixVectorOperations.h"
#include "tarch/la/TransposedMatrix.h"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "utils/Dimensions.hpp"
#include "tarch/la/Scalar.h"
#include "io/TXTWriter.hpp"
#include "io/TXTReader.hpp"
#include "utils/MasterSlave.hpp"
#include "utils/EventTimings.hpp"
#include "utils/EigenHelperFunctions.hpp"
#include "QRFactorization.hpp"
#include "Eigen/Dense"
#include <sys/unistd.h>

#include "tarch/tests/TestMacros.h"

#include <time.h>

//#include "utils/NumericalCompare.hpp"

using precice::utils::Event;

namespace precice {
namespace cplscheme {
namespace impl {

// tarch::logging::Log IQNILSPostProcessing::
//       _log("precice::cplscheme::impl::IQNILSPostProcessing");

IQNILSPostProcessing:: IQNILSPostProcessing
(
  double initialRelaxation,
  bool forceInitialRelaxation,
  int    maxIterationsUsed,
  int    timestepsReused,
  int 	 filter,
  double singularityLimit,
  std::vector<int> dataIDs,
  PtrPreconditioner preconditioner)
:
  BaseQNPostProcessing(initialRelaxation, forceInitialRelaxation, maxIterationsUsed, timestepsReused,
		       filter, singularityLimit, dataIDs, preconditioner),
  _secondaryOldXTildes(),
  _secondaryMatricesW(),
  _secondaryMatricesWBackup()
{
}

void IQNILSPostProcessing:: initialize
(
  DataMap& cplData )
{
  // do common QN post processing initialization
  BaseQNPostProcessing::initialize(cplData);

  double init = 0.0;
  // Fetch secondary data IDs, to be relaxed with same coefficients from IQN-ILS
  for (DataMap::value_type& pair: cplData){
	if (not utils::contained(pair.first, _dataIDs)){
	  int secondaryEntries = pair.second->values->size();
	  utils::append(_secondaryOldXTildes[pair.first], Eigen::VectorXd::Zeros(secondaryEntries));
	}
  }
}


void IQNILSPostProcessing::updateDifferenceMatrices
(
  DataMap& cplData)
{
  Event e(__func__, true, true); // time measurement, barrier
	// Compute residuals of secondary data
	for (int id: _secondaryDataIDs){
		Eigen::VectorXd& secResiduals = _secondaryResiduals[id];
		PtrCouplingData data = cplData[id];
		assertion2(secResiduals.size() == data->values->size(),
				secResiduals.size(), data->values->size());
		secResiduals = *(data->values);
		secResiduals -= data->oldValues.col(0);
	}

	if (_firstIteration && (_firstTimeStep || _forceInitialRelaxation)){
		// constant relaxation: for secondary data called from base class
	}else{
		if (not _firstIteration) {
			bool columnLimitReached = getLSSystemCols() == _maxIterationsUsed;
			bool overdetermined = getLSSystemCols() <= getLSSystemRows();
			if (not columnLimitReached && overdetermined) {

				// Append column for secondary W matrices
				for (int id: _secondaryDataIDs) {
				  utils::appendFront(_secondaryMatricesW[id], _secondaryResiduals[id]);
				}
			}
			else {
				// Shift column for secondary W matrices
				for (int id: _secondaryDataIDs) {
				  utils::shiftSetFirst(_secondaryMatricesW[id], _secondaryResiduals[id]);
				}
			}

			// Compute delta_x_tilde for secondary data
			for (int id: _secondaryDataIDs) {
				Eigen::MatrixXd& secW = _secondaryMatricesW[id];
				assertion2(secW.rows() == cplData[id]->values->size(), secW.rows(), cplData[id]->values->size());
				secW.col(0) = *(cplData[id]->values);
				secW.col(0) -= _secondaryOldXTildes[id];
			}
		}

		// Store x_tildes for secondary data
		for (int id: _secondaryDataIDs) {
			assertion2(_secondaryOldXTildes[id].size() == cplData[id]->values->size(),
					_secondaryOldXTildes[id].size(), cplData[id]->values->size());
			_secondaryOldXTildes[id] = *(cplData[id]->values);
		}
	}
  
  // call the base method for common update of V, W matrices
  BaseQNPostProcessing::updateDifferenceMatrices(cplData);
}


void IQNILSPostProcessing::computeUnderrelaxationSecondaryData
(
  DataMap& cplData)
{
    //Store x_tildes for secondary data
    for (int id: _secondaryDataIDs){
      assertion2(_secondaryOldXTildes[id].size() == cplData[id]->values->size(),
                 _secondaryOldXTildes[id].size(), cplData[id]->values->size());
      _secondaryOldXTildes[id] = *(cplData[id]->values);
    }

    // Perform underrelaxation with initial relaxation factor for secondary data
    for (int id: _secondaryDataIDs){
      PtrCouplingData data = cplData[id];
      Eigen::VectorXd& values = *(data->values);
      values *= _initialRelaxation;                   // new * omg
      Eigen::VectorXd& secResiduals = _secondaryResiduals[id];
      secResiduals = data->oldValues.col(0);    // old
      secResiduals *= 1.0 - _initialRelaxation;       // (1-omg) * old
      values += secResiduals;                      // (1-omg) * old + new * omg
    }
}


void IQNILSPostProcessing::computeQNUpdate
(PostProcessing::DataMap& cplData, Eigen::VectorXd& xUpdate)
{
	preciceTrace("computeQNUpdate()");
  Event e(__func__, true, true); // time measurement, barrier
  preciceDebug("   Compute Newton factors");

  // Calculate QR decomposition of matrix V and solve Rc = -Qr
  Eigen::VectorXd __c;

	// for master-slave mode and procs with no vertices,
	// qrV.cols() = getLSSystemCols() and _qrV.rows() = 0
	auto __Qt = _qrV.matrixQ();
	auto __R = _qrV.matrixR();

	if(!_hasNodesOnInterface){
	  assertion2(_qrV.cols() == getLSSystemCols(), _qrV.cols(), getLSSystemCols());
	  assertion1(_qrV.rows() == 0, _qrV.rows());
	  assertion1(__Qt.size() == 0, __Qt.size());
	}

	Eigen::VectorXd _local_b(_qrV.cols(), 0.0);
	Eigen::VectorXd _global_b;

	Event e_qrsolve("solve: R alpha = -Q^T r", true, true); // time measurement, barrier
	multiply(__Qt, _residuals, _local_b);
	_local_b *= -1.0; // = -Qr

	assertion1(__c.size() == 0, __c.size());
	// reserve memory for c
	utils::append(__c, Eigen::VectorXd::Zero(_local_b.size()));

	// compute rhs Q^T*res in parallel
	if (not utils::MasterSlave::_masterMode && not utils::MasterSlave::_slaveMode) {
		assertion2(__Qt.rows() == getLSSystemCols(), __Qt.rows(), getLSSystemCols());
		// back substitution
		__c = __R.triangularView<Eigen::Upper>().solve<Eigen::OnTheLeft>(_local_b);
	}else{
	   assertion(utils::MasterSlave::_communication.get() != nullptr);
	   assertion(utils::MasterSlave::_communication->isConnected());
	   if(_hasNodesOnInterface)  assertion2(__Qt.rows() == getLSSystemCols(), __Qt.rows(), getLSSystemCols());
	   assertion2(_local_b.size() == getLSSystemCols(), _local_b.size(), getLSSystemCols());

	   if(utils::MasterSlave::_masterMode){
	     assertion1(_global_b.size() == 0, _global_b.size());
	   }
	   utils::append(_global_b, Eigen::VectorXd::Zero(_local_b.size()));

	   // do a reduce operation to sum up all the _local_b vectors
	   utils::MasterSlave::reduceSum(_local_b.data(), _global_b.data(), _local_b.size()); // size = getLSSystemCols() = _local_b.size()

	   // back substitution R*c = b only in master node
	   if(utils::MasterSlave::_masterMode)
	     __c = __R.triangularView<Eigen::Upper>().solve<Eigen::OnTheLeft>(_global_b);

	  // broadcast coefficients c to all slaves
	  utils::MasterSlave::broadcast(__c.data(), __c.size());
	}
	e_qrsolve.stop();

	preciceDebug("   Apply Newton factors");
	// compute x updates from W and coefficients c, i.e, xUpdate = c*W
	xUpdate = _matrixW * __c;

	preciceDebug("c = " << __c);


    /**
     *  perform QN-Update step for the secondary Data
     */

	// If the previous time step converged within one single iteration, nothing was added
	// to the LS system matrices and they need to be restored from the backup at time T-2
    if (not _firstTimeStep && (getLSSystemCols() < 1) && (_timestepsReused == 0) && not _forceInitialRelaxation) {
		preciceDebug("   Last time step converged after one iteration. Need to restore the secondaryMatricesW from backup.");
		_secondaryMatricesW = _secondaryMatricesWBackup;
	}

	// Perform QN relaxation for secondary data
	for (int id: _secondaryDataIDs){
	  PtrCouplingData data = cplData[id];
	  Eigen::VectorXd& values = *(data->values);
	  assertion2(_secondaryMatricesW[id].cols() == __c.size(), _secondaryMatricesW[id].cols(), __c.size());
	  values = _secondaryMatricesW[id] * __c;
	  assertion2(values.size() == data->oldValues.column(0).size(), values.size(), data->oldValues.column(0).size());
	  values += data->oldValues.col(0);
	  assertion2(values.size() == _secondaryResiduals[id].size(), values.size(), _secondaryResiduals[id].size());
	  values += _secondaryResiduals[id];
	}

	// pending deletion: delete old secondaryMatricesW
	if (_firstIteration && _timestepsReused == 0 && not _forceInitialRelaxation) {
		// save current secondaryMatrix data in case the coupling for the next time step will terminate
		// after the first iteration (no new data, i.e., V = W = 0)
		if(getLSSystemCols() > 0){
			_secondaryMatricesWBackup = _secondaryMatricesW;
		}
		for (int id: _secondaryDataIDs){
			_secondaryMatricesW[id].resize(0,0);
		}
	}
}



void IQNILSPostProcessing:: specializedIterationsConverged
(
   DataMap & cplData)
{
  
  if (_matrixCols.front() == 0){ // Did only one iteration
    _matrixCols.pop_front(); 
  }

  if (_timestepsReused == 0){
    if (_forceInitialRelaxation)
    {
      for (int id: _secondaryDataIDs) {
        _secondaryMatricesW[id].resize(0,0);
      }
    } else {
      /**
       * pending deletion (after first iteration of next time step
       * Using the matrices from the old time step for the first iteration
       * is better than doing underrelaxation as first iteration of every time step
       */
    }
  }
  else if ((int)_matrixCols.size() > _timestepsReused){
    int toRemove = _matrixCols.back();
    for (int id: _secondaryDataIDs){
      Eigen::MatrixXd& secW = _secondaryMatricesW[id];
      assertion3(secW.cols() > toRemove, secW, toRemove, id);
      for (int i=0; i < toRemove; i++){
        utils::removeColumnFromMatrix(secW, secW.cols() - 1);
      }
    }
  }
}


void IQNILSPostProcessing:: removeMatrixColumn
(
  int columnIndex)
{
  assertion(_matrixV.cols() > 1);
  // remove column from secondary Data Matrix W
  for (int id: _secondaryDataIDs){
    utils::removeColumnFromMatrix(_secondaryMatricesW[id], columnIndex);
   }

	BaseQNPostProcessing::removeMatrixColumn(columnIndex);
}

}}} // namespace precice, cplscheme, impl
