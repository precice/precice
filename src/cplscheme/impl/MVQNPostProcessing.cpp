// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "MVQNPostProcessing.hpp"
#include "cplscheme/CouplingData.hpp"
#include "utils/Globals.hpp"
#include "tarch/la/GramSchmidt.h"
#include "tarch/la/LUDecomposition.h"
#include "tarch/la/MatrixVectorOperations.h"
#include "tarch/la/TransposedMatrix.h"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "utils/Dimensions.hpp"
#include "utils/MasterSlave.hpp"
#include "tarch/la/Scalar.h"
#include "io/TXTWriter.hpp"
#include "io/TXTReader.hpp"
#include "Eigen/Dense"

#include <time.h>
#include <sstream>
#include <fstream>
//#include "utils/NumericalCompare.hpp"

namespace precice {
namespace cplscheme {
namespace impl {

// tarch::logging::Log MVQNPostProcessing::
//       _log("precice::cplscheme::impl::MVQNPostProcessing");

      
MVQNPostProcessing:: MVQNPostProcessing
(
  double initialRelaxation,
  int    maxIterationsUsed,
  int    timestepsReused,
  double singularityLimit,
  std::vector<int> dataIDs,
  std::map<int,double> scalings)
:
  BaseQNPostProcessing(initialRelaxation, maxIterationsUsed, timestepsReused,
		       singularityLimit, dataIDs, scalings),
//  _secondaryOldXTildes(),
  _invJacobian(),
  _oldInvJacobian(),
  _dimOffsets()
{}



void MVQNPostProcessing:: initialize
(
  DataMap& cplData )
{
  // do common QN post processing initialization
  BaseQNPostProcessing::initialize(cplData);
  
  double init = 0.0;
  size_t entries = _residuals.size();
  size_t global_n = 0;

	if (not utils::MasterSlave::_masterMode && not utils::MasterSlave::_slaveMode) {
		global_n = entries;
	}else{

		assertion(utils::MasterSlave::_communication.get() != NULL);
		assertion(utils::MasterSlave::_communication->isConnected());

		/**
		 *  make dimensions public to all procs,
		 *  last entry _dimOffsets[MasterSlave::_size] holds the global dimension, global,n
		 */
		_dimOffsets.resize(utils::MasterSlave::_size + 1);
		if (utils::MasterSlave::_slaveMode) {
			utils::MasterSlave::_communication->send(_matrixV.rows(), 0);
			utils::MasterSlave::_communication->receive(&_dimOffsets[0], _dimOffsets.size(), 0);
		}
		if (utils::MasterSlave::_masterMode) {
			_dimOffsets[0] = 0;
			_dimOffsets[1] = _matrixV.rows();
			for (int rankSlave = 1; rankSlave < _size; rankSlave++) {
				int localDim = 0;
				utils::MasterSlave::_communication->receive(localDim, rankSlave);
				_dimOffsets[rankSlave + 1] = _dimOffsets[rankSlave] + localDim;
			}
			for (int rankSlave = 1; rankSlave < _size; rankSlave++) {
				utils::MasterSlave::_communication->send(&_dimOffsets[0], _dimOffsets.size(), rankSlave);
			}
		}
		global_n = _dimOffsets.back();
	}
  
  _invJacobian = Matrix(global_n, entries, init);
  _oldInvJacobian = Matrix(global_n, entries, init);
}



void MVQNPostProcessing::computeUnderrelaxationSecondaryData
(
  DataMap& cplData)
{
    //Store x_tildes for secondary data
  //  foreach (int id, _secondaryDataIDs){
  //    assertion2(_secondaryOldXTildes[id].size() == cplData[id]->values->size(),
  //               _secondaryOldXTildes[id].size(), cplData[id]->values->size());
  //    _secondaryOldXTildes[id] = *(cplData[id]->values);
  //  }

    // Perform underrelaxation with initial relaxation factor for secondary data
    foreach (int id, _secondaryDataIDs){
      PtrCouplingData data = cplData[id];
      DataValues& values = *(data->values);
      values *= _initialRelaxation;                   // new * omg
      DataValues& secResiduals = _secondaryResiduals[id];
      secResiduals = data->oldValues.column(0);    // old
      secResiduals *= 1.0 - _initialRelaxation;       // (1-omg) * old
      values += secResiduals;                      // (1-omg) * old + new * omg
    }
}




void MVQNPostProcessing::updateDifferenceMatrices
(
  DataMap& cplData)
{
  using namespace tarch::la;

//   // Compute residuals of secondary data
//   foreach (int id, _secondaryDataIDs){
//     DataValues& secResiduals = _secondaryResiduals[id];
//     PtrCouplingData data = cplData[id];
//     assertion2(secResiduals.size() == data->values->size(),
//                secResiduals.size(), data->values->size());
//     secResiduals = *(data->values);
//     secResiduals -= data->oldValues.column(0);
//   }

  /*
   * ATTETION: changed the condition from _firstIteration && _firstTimeStep
   * to the following: 
   * underrelaxation has to be done, if the scheme has converged without even
   * entering post processing. In this case the V, W matrices would still be empty.
   * This case happended in the open foam example beamInCrossFlow.
   */ 
  if(_firstIteration && (_firstTimeStep ||  (_matrixCols.size() < 2))){
    //k++;
    // Perform underrelaxation with initial relaxation factor for secondary data
//     foreach (int id, _secondaryDataIDs){
//       PtrCouplingData data = cplData[id];
//       DataValues& values = *(data->values);
//       values *= _initialRelaxation;                   // new * omg
//       DataValues& secResiduals = _secondaryResiduals[id];
//       secResiduals = data->oldValues.column(0);    // old
//       secResiduals *= 1.0 - _initialRelaxation;       // (1-omg) * old
//       values += secResiduals;                      // (1-omg) * old + new * omg
//     }
  }
  else {
    if (not _firstIteration){
      //k++;
    }
  }
  
  // call the base method for common update of V, W matrices
  BaseQNPostProcessing::updateDifferenceMatrices(cplData);
}



void MVQNPostProcessing::computeQNUpdate
    (PostProcessing::DataMap& cplData, DataValues& xUpdate)
{
  preciceTrace("computeQNUpdate()");
  using namespace tarch::la;

    // ------------- update inverse Jacobian -----------
    // J_inv = J_inv_n + (W - J_inv_n*V)*(V^T*V)^-1*V^T
    // ----------------------------------------- -------

    preciceDebug("   Compute Newton factors");
    //computeNewtonFactorsLUDecomposition(cplData, xUpdate);
    
    // computes xUpdate using updatedQR decompositon (does not modify _invJacobian)
    computeNewtonFactorsUpdatedQRDecomposition(cplData, xUpdate);
}


void MVQNPostProcessing::computeNewtonFactorsUpdatedQRDecomposition
(PostProcessing::DataMap& cplData, DataValues& xUpdate)
{
  preciceTrace("computeNewtonFactorsQRDecomposition()");
  using namespace tarch::la;
 
  // ------------- update inverse Jacobian -----------
  // J_inv = J_inv_n + (W - J_inv_n*V)*(V^T*V)^-1*V^T
  // ----------------------------------------- -------

  DataMatrix Z;
  bool linearDependence = true;
  
  while (linearDependence) {
		linearDependence = false;
		Z.clear();

		Matrix __R(_matrixV.cols(), _matrixV.cols(), 0.0);
		auto r = _qrV.matrixR();
		for (int i = 0; i < r.rows(); i++)
			for (int j = 0; j < r.cols(); j++) {
				__R(i, j) = r(i, j);
			}

		if (_matrixV.cols() > 1) {
			for (int i = 0; i < _matrixV.cols(); i++) {
				if (std::fabs(__R(i, i)) < _singularityLimit) {
					preciceDebug("   Removing linear dependent column " << i);
					_infostream
							<< "(updatedQR) removing linear dependent column "
							<< i << "  time step: " << tSteps << " iteration: " << its
							<< "\n" << std::flush;
					linearDependence = true;
					removeMatrixColumn(i);
				}
			}
		}
		if (not linearDependence) {
			Matrix __Q(_matrixV.rows(), _matrixV.cols(), 0.0);

			DataValues __ytmpVec(_matrixV.cols(), 0.0);
			DataValues __matrixQRow;
			auto q = _qrV.matrixQ();
			for (int i = 0; i < q.rows(); i++)
				for (int j = 0; j < q.cols(); j++) {
					__Q(i, j) = q(i, j);
				}

			r = _qrV.matrixR();
			for (int i = 0; i < r.rows(); i++)
				for (int j = 0; j < r.cols(); j++) {
					__R(i, j) = r(i, j);
				}
			for (int i = 0; i < __Q.rows(); i++) {
				for (int j = 0; j < __Q.cols(); j++) {
					__matrixQRow.append(__Q(i, j));
				}

				backSubstitution(__R, __matrixQRow, __ytmpVec);
				Z.append(__ytmpVec);
				__matrixQRow.clear();
			}
		}
	}


  /*
   * Multiply J_prev * V =: V_tilde
   */
  // TODO: transpose V efficiently using blocking in parallel
  //       such that multiplication is cache efficient
  Matrix tmpMatrix(_matrixV.rows(), _matrixV.cols(), 0.0);
  assertion2(_oldInvJacobian.cols() == _matrixV.rows(), _oldInvJacobian.cols(), _matrixV.rows());

  if (not utils::MasterSlave::_masterMode && not utils::MasterSlave::_slaveMode) {
	  multiply(_oldInvJacobian, _matrixV, tmpMatrix);

  }else{

	  assertion(utils::MasterSlave::_communication.get() != NULL);
	  assertion(utils::MasterSlave::_communication->isConnected());

	  for(int i = 0; i < tmpMatrix.rows(); i++){
		  // find row of _oldInvJacobian that corresponds to processor block
		  int jrow = offsets[utils::MasterSlave::_rank] + i;
		  for(int j = 0; j < tmpMatrix.cols(); j++){
			  // as we want to move to Eigen, copy
			  DataValues Jrow(_oldInvJacobian.cols(), 0.0);
			  for (int s = 0; s < _oldInvJacobian.cols(); s++) {
				  Jrow(s) = _oldInvJacobian(jrow,s);
			  }
			  // TODO: better: implement a reduce-operation (no loop over all slaves)
			  tmpMatrix(i,j) = utils::MasterSlave::dot(Jrow, _matrixV.column(j));
		  }
	  }
  }

  // tmpMatrix = (W-J_inv_n*V)
  tmpMatrix *= -1.;
  tmpMatrix = tmpMatrix + _matrixW;


  
  /**
   *  compute invJacobian = W_til*Z
   *  where Z = (V^T*V)^-1*V^T vie QR-dec and back-substitution
   *  and W_til = (W - J_inv_n*V)
   */
  assertion2(tmpMatrix.cols() == Z.rows(), tmpMatrix.cols(), Z.rows());
  if (not utils::MasterSlave::_masterMode && not utils::MasterSlave::_slaveMode) {
	  multiply(tmpMatrix, Z, _invJacobian);

  }else{

	  assertion(utils::MasterSlave::_communication.get() != NULL);
	  assertion(utils::MasterSlave::_communication->isConnected());

	  // compute block matrices of J_inv of size (n_til x n_til), n_til = local n
	  Matrix diagBlock(_matrixV.rows(),_matrixV.rows(), 0.0);
	  multiply(tmpMatrix, Z, diagBlock);
	  // set block at corresponding row-index on proc
	  int off = offsets[utils::MasterSlave::_rank];
	  assertion2(_invJacobian.cols() == diagBlock.cols(), _invJacobian.cols(), diagBlock.cols());
	  for(int q = 0; q < diagBlock.rows(); q++)
		  for(int p = 0; p < _invJacobian.cols(); p++)
		  {
			  _invJacobian(q+off,p) = diagBlock(q,p);
		  }

	  for(int cycle = 1; cycle < utils::MasterSlave::_size; cycle++){
		  // cyclic send-receive operation
		  int nextProc = (utils::MasterSlave::_rank + 1) % utils::MasterSlave::_size;
		  int prevProc = (utils::MasterSlave::_rank -1 < 0) ? utils::MasterSlave::_size-1 : utils::MasterSlave::_rank -1;
		  int rows_rcv = (prevProc > 0) ? offsets[prevProc] - offsets[prevProc-1] : offsets[0];
		  Matrix Wtil_rcv(rows_rcv, tmpMatrix.cols(),0.0);

		  utils::MasterSlave::_communication->send(&tmpMatrix(0,0), tmpMatrix.rows(), tmpMatrix.cols(), nextProc); //TODO: implementation
		  utils::MasterSlave::_communication->receive(&Wtil_rcv(0,0), rows_rcv, tmpMatrix.cols(), prevProc);

		  // compute block with new local data
		  Matrix block(rows_rcv, Z.cols(), 0.0);
		  multiply(Wtil_rcv, Z, block);
		  // set block at corresponding index in J_inv
		  // the row-offset of the current block is determined by the proc that sends the part of the W_til matrix
		  // note: the direction and ordering of the cyclic sending operation is chosen s.t. the computed block is
		  //       local on the current processor (in J_inv).
		  off = offsets[prevProc];
		  assertion2(_invJacobian.cols() == block.cols(), _invJacobian.cols(), block.cols());
		  for(int q = 0; q < block.rows(); q++)
			  for(int p = 0; p < _invJacobian.cols(); p++)
			  {
				  _invJacobian(q+off,p) = block(q,p);
			  }
	  }
  }

  _invJacobian = _invJacobian + _oldInvJacobian;

  DataValues negRes(_residuals);
  negRes *= -1.;
  /*
   * solve delta_x = - J_inv*residuals
   *
   * two variants possible:
   * (1) Similar to the computation of V_til = J_inv_n*V loop over rows
   *     and cols of xUpdates = J_inv*res and compute MasterSlave::dot(J.row(i), res)
   * (2) Similar to computation of rhs in IQN-ILS multiply local block J_inv*res,
   *     save result in tmp vector of size n_global and send it to the master.
   *
   *  First variant is maybe better, as for the second variant, each tmp result
   *  has to be published to all slaves such that they can compute their part of
   *  th xUpdate
   */
  if (not utils::MasterSlave::_masterMode && not utils::MasterSlave::_slaveMode) {
	  multiply(_invJacobian, negRes, xUpdate);

  }else{

  	  assertion(utils::MasterSlave::_communication.get() != NULL);
  	  assertion(utils::MasterSlave::_communication->isConnected());

  	  for(int i = 0; i < xUpdate.size(); i++){
  		  // find row of _oldInvJacobian that corresponds to processor block
  		  int jrow = offsets[utils::MasterSlave::_rank] + i;
		  // as we want to move to Eigen, copy
		  DataValues Jrow(_invJacobian.cols(), 0.0);
		  for (int s = 0; s < _invJacobian.cols(); s++) {
			  Jrow(s) = _invJacobian(jrow,s);
		  }
		  // TODO: better: implement a reduce-operation (no loop over all slaves)
		  xUpdate(i) = utils::MasterSlave::dot(Jrow, negRes);
  	  }
    }
}


void MVQNPostProcessing::computeNewtonFactorsQRDecomposition
(PostProcessing::DataMap& cplData, DataValues& xUpdate)
{
  preciceTrace("computeNewtonFactorsQRDecomposition()");
  using namespace tarch::la;
 
  // ------------- update inverse Jacobian -----------
  // J_inv = J_inv_n + (W - J_inv_n*V)*(V^T*V)^-1*V^T
  // ----------------------------------------- -------
  
  DataMatrix v;
  bool linearDependence = true;
	while (linearDependence) {
		linearDependence = false;
		v.clear();

		DataMatrix Vcopy(_matrixV);
		DataMatrix Q(Vcopy.rows(), Vcopy.cols(), 0.0);
		DataMatrix R(Vcopy.cols(), Vcopy.cols(), 0.0);

		modifiedGramSchmidt(Vcopy, Q, R);

		if (_matrixV.cols() > 1) {
			for (int i = 0; i < _matrixV.cols(); i++) {
				if (std::fabs(R(i, i)) < _singularityLimit) {
					preciceDebug("   Removing linear dependent column " << i);
					_infostream
							<< "(modifiedGramSchmidt) removing linear dependent column "
							<< i << "  time step: " << tSteps << " iteration: " << its
							<< "\n" << std::flush;
					linearDependence = true;
					removeMatrixColumn(i);
				}
			}
		}
		if (not linearDependence) {
			DataValues ytmpVec(_matrixV.cols(), 0.0);
			DataValues _matrixQRow;
			for (int i = 0; i < Q.rows(); i++) {
				for (int j = 0; j < Q.cols(); j++) {
					_matrixQRow.append(Q(i, j));
				}
				backSubstitution(R, _matrixQRow, ytmpVec);
				v.append(ytmpVec);
				_matrixQRow.clear();
			}
		}
	}

  // tmpMatrix = J_inv_n*V
  Matrix tmpMatrix(_matrixV.rows(), _matrixV.cols(), 0.0);
  assertion2(_oldInvJacobian.cols() == _matrixV.rows(), _oldInvJacobian.cols(), _matrixV.rows());
  multiply(_oldInvJacobian, _matrixV, tmpMatrix);

  // tmpMatrix = (W-J_inv_n*V)
  tmpMatrix *= -1.;
  tmpMatrix = tmpMatrix + _matrixW;
  
  // invJacobian = (W - J_inv_n*V)*(V^T*V)^-1*V^T
  assertion2(tmpMatrix.cols() == v.rows(), tmpMatrix.cols(), v.rows());
  Matrix tmp_invJacobian(_invJacobian.rows(), _invJacobian.cols(), 0.0);
  multiply(tmpMatrix, v, tmp_invJacobian);
  tmp_invJacobian = tmp_invJacobian + _oldInvJacobian;
  
  DataValues negRes(_residuals);
  negRes *= -1.;
  
  // solve delta_x = - J_inv*residuals
  multiply(tmp_invJacobian, negRes, xUpdate); 
}


void MVQNPostProcessing::computeNewtonFactorsLUDecomposition
(PostProcessing::DataMap& cplData, DataValues& xUpdate)
{
  preciceTrace("computeNewtonFactorsLUDecomposition()");
  using namespace tarch::la;
  
  // ------------- update inverse Jacobian -----------
  // J_inv = J_inv_n + (W - J_inv_n*V)*(V^T*V)^-1*V^T
  // ----------------------------------------- -------
  
  DataMatrix VTVLU(_matrixV.cols(), _matrixV.cols(), 0.0);
  DataMatrix v;
  multiply(transpose(_matrixV), _matrixV, VTVLU);  // VTV = V^T*
  
  DataValues pivots(_matrixV.cols(), 0.0);
  lu(VTVLU,pivots);
  
  
  DataValues ytmpVec(_matrixV.cols(), 0.0);
  DataValues xtmpVec(_matrixV.cols(), 0.0);
  DataValues _matrixVRow;
  for(int i = 0; i < _matrixV.rows(); i++)
  {
    for(int j=0; j < _matrixV.cols(); j++){
      _matrixVRow.append(_matrixV(i,j));
    }
    
    // account for pivoting in lu-decomposition
    assertion2(_matrixVRow.size() == pivots.size(), _matrixVRow.size(), pivots.size());
    for ( int i=0; i < _matrixVRow.size(); i++ ){
      double temp = _matrixVRow[i];
      _matrixVRow[i] = _matrixVRow[pivots[i]];
      _matrixVRow[pivots[i]] = temp;
    }
    forwardSubstitution(VTVLU, _matrixVRow, ytmpVec);

    backSubstitution(VTVLU, ytmpVec, xtmpVec);
    
    v.append(xtmpVec);  
    _matrixVRow.clear();
  }
  
  // tmpMatrix = J_inv_n*V
  DataMatrix tmpMatrix(_matrixV.rows(), _matrixV.cols(), 0.0);
  assertion2(_oldInvJacobian.cols() == _matrixV.rows(), _oldInvJacobian.cols(), _matrixV.rows());
  multiply(_oldInvJacobian, _matrixV, tmpMatrix);
  // tmpMatrix = (W-J_inv_n*V)
  tmpMatrix *= -1.;
  tmpMatrix = _matrixW + tmpMatrix;
  // invJacobian = (W - J_inv_n*V)*(V^T*V)^-1*V^T
  assertion2(tmpMatrix.cols() == v.rows(), tmpMatrix.cols(), v.rows());
  multiply(tmpMatrix, v, _invJacobian);
  _invJacobian = _invJacobian + _oldInvJacobian;
  
  DataValues negRes(_residuals);
  negRes *= -1.;
  
  // solve delta_x = - J_inv*residuals
  multiply(_invJacobian, negRes, xUpdate); 
}



void MVQNPostProcessing:: specializedIterationsConverged
(
   DataMap & cplData)
{
  
  
  //k = 0;
  //t++;
  // store inverse Jacobian
//  _matrixWriter.write(_invJacobian);
  _oldInvJacobian = _invJacobian;
}

}}} // namespace precice, cplscheme, impl
