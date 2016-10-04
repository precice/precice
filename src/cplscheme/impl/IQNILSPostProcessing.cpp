#include "IQNILSPostProcessing.hpp"
#include "cplscheme/CouplingData.hpp"
#include "utils/Globals.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "utils/Dimensions.hpp"
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

// logging::Logger IQNILSPostProcessing::
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

void IQNILSPostProcessing::initialize
(
  DataMap& cplData )
{
  Event e("initialize()", true, true); // time measurement, barrier

  // do common QN post processing initialization
  BaseQNPostProcessing::initialize(cplData);

  double init = 0.0;
  // Fetch secondary data IDs, to be relaxed with same coefficients from IQN-ILS
  for (DataMap::value_type& pair: cplData){
    if (not utils::contained(pair.first, _dataIDs)){
      int secondaryEntries = pair.second->values->size();
      utils::append(_secondaryOldXTildes[pair.first], (Eigen::VectorXd) Eigen::VectorXd::Zero(secondaryEntries));
    }
  }
}


void IQNILSPostProcessing::updateDifferenceMatrices
(
  DataMap& cplData)
{
  Event e("IQNILSPostProcessing::updateDifferenceMatrices", true, true); // time measurement, barrier
	// Compute residuals of secondary data
	for (int id: _secondaryDataIDs){
		Eigen::VectorXd& secResiduals = _secondaryResiduals[id];
		PtrCouplingData data = cplData[id];
		assertion(secResiduals.size() == data->values->size(),
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
				assertion(secW.rows() == cplData[id]->values->size(), secW.rows(), cplData[id]->values->size());
				secW.col(0) = *(cplData[id]->values);
				secW.col(0) -= _secondaryOldXTildes[id];
			}
		}

		// Store x_tildes for secondary data
		for (int id: _secondaryDataIDs) {
			assertion(_secondaryOldXTildes[id].size() == cplData[id]->values->size(),
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
      assertion(_secondaryOldXTildes[id].size() == cplData[id]->values->size(),
                 _secondaryOldXTildes[id].size(), cplData[id]->values->size());
      _secondaryOldXTildes[id] = *(cplData[id]->values);
    }

    // Perform underrelaxation with initial relaxation factor for secondary data
    for (int id: _secondaryDataIDs){
      PtrCouplingData data = cplData[id];
      auto& values = *(data->values);
      auto& secResiduals = _secondaryResiduals[id];
      secResiduals = (1.0 - _initialRelaxation) * data->oldValues.col(0);    // old

      values +=  (values * _initialRelaxation) + secResiduals;
    }
}


void IQNILSPostProcessing::computeQNUpdate
(PostProcessing::DataMap& cplData, Eigen::VectorXd& xUpdate)
{
	preciceTrace("computeQNUpdate()");
  Event e("computeNewtonUpdate", true, true); // time measurement, barrier

  DEBUG("   Compute Newton factors");

  // Calculate QR decomposition of matrix V and solve Rc = -Qr
  Eigen::VectorXd c;

	// for master-slave mode and procs with no vertices,
	// qrV.cols() = getLSSystemCols() and _qrV.rows() = 0
	auto Q = _qrV.matrixQ();
	auto R = _qrV.matrixR();

	if(!_hasNodesOnInterface){
	  assertion(_qrV.cols() == getLSSystemCols(), _qrV.cols(), getLSSystemCols());
	  assertion(_qrV.rows() == 0, _qrV.rows());
	  assertion(Q.size() == 0, Q.size());
	}

	Eigen::VectorXd _local_b = Eigen::VectorXd::Zero(_qrV.cols());
	Eigen::VectorXd _global_b;

	Event e_qrsolve("solve: R alpha = -Q^T r", true, true); // time measurement, barrier

	// need to scale the residual to compensate for the scaling in c = R^-1 * Q^T * P^-1 * residual'
	// it is also possible to apply the inverse scaling weights from the right to the vector c
	_preconditioner->apply(_residuals);
	_local_b = Q.transpose() * _residuals;
	_preconditioner->revert(_residuals);
	_local_b *= -1.0; // = -Qr

	assertion(c.size() == 0, c.size());
	// reserve memory for c
	utils::append(c, (Eigen::VectorXd) Eigen::VectorXd::Zero(_local_b.size()));

	// compute rhs Q^T*res in parallel
	if (not utils::MasterSlave::_masterMode && not utils::MasterSlave::_slaveMode) {
		assertion(Q.cols() == getLSSystemCols(), Q.cols(), getLSSystemCols());
		// back substitution
		c = R.triangularView<Eigen::Upper>().solve<Eigen::OnTheLeft>(_local_b);
	}else{
	   assertion(utils::MasterSlave::_communication.get() != nullptr);
	   assertion(utils::MasterSlave::_communication->isConnected());
	   if(_hasNodesOnInterface)  assertion(Q.cols() == getLSSystemCols(), Q.cols(), getLSSystemCols());
	   assertion(_local_b.size() == getLSSystemCols(), _local_b.size(), getLSSystemCols());

	   if(utils::MasterSlave::_masterMode){
	     assertion(_global_b.size() == 0, _global_b.size());
	   }
	   utils::append(_global_b, (Eigen::VectorXd) Eigen::VectorXd::Zero(_local_b.size()));

	   // do a reduce operation to sum up all the _local_b vectors
	   utils::MasterSlave::reduceSum(_local_b.data(), _global_b.data(), _local_b.size()); // size = getLSSystemCols() = _local_b.size()

	   // back substitution R*c = b only in master node
	   if(utils::MasterSlave::_masterMode)
	     c = R.triangularView<Eigen::Upper>().solve<Eigen::OnTheLeft>(_global_b);

	  // broadcast coefficients c to all slaves
	  utils::MasterSlave::broadcast(c.data(), c.size());
	}
	e_qrsolve.stop();

	// REOMVE!!
	//c = _matrixV.householderQr().solve(-_residuals);


	DEBUG("   Apply Newton factors");
	// compute x updates from W and coefficients c, i.e, xUpdate = c*W
	//xUpdate = ((_matrixW-_matrixV) + 1.*_matrixV) * c;
	xUpdate = _matrixW * c;

	//DEBUG("c = " << c);


  /**
   *  perform QN-Update step for the secondary Data
   */

	// If the previous time step converged within one single iteration, nothing was added
	// to the LS system matrices and they need to be restored from the backup at time T-2
    if (not _firstTimeStep && (getLSSystemCols() < 1) && (_timestepsReused == 0) && not _forceInitialRelaxation) {
		DEBUG("   Last time step converged after one iteration. Need to restore the secondaryMatricesW from backup.");
		_secondaryMatricesW = _secondaryMatricesWBackup;
	}

	// Perform QN relaxation for secondary data
	for (int id: _secondaryDataIDs){
	  PtrCouplingData data = cplData[id];
	  auto& values = *(data->values);
	  assertion(_secondaryMatricesW[id].cols() == c.size(), _secondaryMatricesW[id].cols(), c.size());
	  values = _secondaryMatricesW[id] * c;
	  assertion(values.size() == data->oldValues.col(0).size(), values.size(), data->oldValues.col(0).size());
	  values += data->oldValues.col(0);
	  assertion(values.size() == _secondaryResiduals[id].size(), values.size(), _secondaryResiduals[id].size());
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
  Event e(__func__, true, true); // time measurement, barrier

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
      assertion(secW.cols() > toRemove, secW, toRemove, id);
      for (int i=0; i < toRemove; i++){
        utils::removeColumnFromMatrix(secW, secW.cols() - 1);
      }
    }
  }
  //e.stop(true);
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
