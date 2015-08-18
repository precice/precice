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
  _oldInvJacobian()
{}



void MVQNPostProcessing:: initialize
(
  DataMap& cplData )
{
  // do common QN post processing initialization
  BaseQNPostProcessing::initialize(cplData);
  
  double init = 0.0;
  int entries = _residuals.size();
  int global_n = 0;

	if (not utils::MasterSlave::_masterMode && not utils::MasterSlave::_slaveMode) {
		global_n = entries;
	}else{
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





/**
 * This is a no-op at the moment, add implementation to handle secondary data
 */
void MVQNPostProcessing::updateDifferenceMatrices
(
  DataMap& cplData)
{
  using namespace tarch::la;

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
    
    // computes xUpdate using updatedQR decompositon
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

  Matrix Z;
  bool linearDependence = true;
  
  while (linearDependence) {
		linearDependence = false;
		//Z.clear();

		Matrix __R(_qrV.cols(), _qrV.cols(), 0.0);
		auto r = _qrV.matrixR();
		for (int i = 0; i < r.rows(); i++)
			for (int j = 0; j < r.cols(); j++) {
				__R(i, j) = r(i, j);
			}

		if(!_hasNodesOnInterface){
						std::cout<<"empty proc["<<utils::MasterSlave::_rank<<"] _qrV: rows="<<_qrV.rows()<<", cols="<<_qrV.cols()<<std::endl;
		}
		if (getLSSystemCols() > 1) {
			if(!_hasNodesOnInterface){
				std::cout<<" proc["<<utils::MasterSlave::_rank<<"] R:\n"<<__R<<std::endl;
			}
			for (int i = 0; i < __R.rows(); i++) {
				//if (std::fabs(__R(i, i)) < _singularityLimit) {
				if (std::fabs(__R(i, i)) < 0.0) {
					std::cout<<" -remove: proc["<<utils::MasterSlave::_rank<<"] R(i,i):"<<__R(i,i)<<", eps:"<<_singularityLimit<<", i:"<<i<<std::endl;
					std::stringstream ss;
					ss << "(updatedQR) removing linear dependent column "
					   << i << "  time step: " << tSteps << " iteration: " << its
					   << "\n" << std::endl;
					preciceDebug(ss.str());
					writeInfo(ss.str());
					std::cout<<ss.str()<<std::endl;

					linearDependence = true;
					removeMatrixColumn(i);
				}
			}
		}
		if (not linearDependence) {
			Matrix __Q(_qrV.rows(), _qrV.cols(), 0.0);
			Z = Matrix(_qrV.cols(), _qrV.rows(), 0.0);

			DataValues __ytmpVec(_qrV.cols(), 0.0);
			DataValues __matrixQRow;
			auto q = _qrV.matrixQ();
			for (int i = 0; i < q.rows(); i++)
				for (int j = 0; j < q.cols(); j++) {
					__Q(i, j) = q(i, j);
				}

			// assertions for the case of processors with no vertices
			if(!_hasNodesOnInterface){
					assertion2(_qrV.cols() == getLSSystemCols(), _qrV.cols(), getLSSystemCols()); assertion1(_qrV.rows() == 0, _qrV.rows()); assertion1(__Q.size() == 0, __Q.size());
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
				for(int p = 0; p < __ytmpVec.size(); p++)
					Z(p,i) = __ytmpVec(p);
				__matrixQRow.clear();
			}
			if(!_hasNodesOnInterface) std::cout<<"\nZ: "<<Z<<std::endl;
			preciceDebug("doing back substitution, __ytmpVec.size="<<__ytmpVec.size());
		}
	}

  /*
   * Multiply J_prev * V =: V_tilde
   */
  // TODO: transpose V efficiently using blocking in parallel
  //       such that multiplication is cache efficient
  Matrix tmpMatrix(_qrV.rows(), _qrV.cols(), 0.0);
  assertion2(_matrixV.rows() == _qrV.rows(), _matrixV.rows(), _qrV.rows());
  assertion2(_oldInvJacobian.cols() == _matrixV.rows(), _oldInvJacobian.cols(), _matrixV.rows());

  preciceDebug("multiply J_prev * V | _invJacobian rows:"<<_oldInvJacobian.rows()<<", cols: "<<_oldInvJacobian.cols());
  if (not utils::MasterSlave::_masterMode && not utils::MasterSlave::_slaveMode) {
	  multiply(_oldInvJacobian, _matrixV, tmpMatrix);

  }else{

	  assertion(utils::MasterSlave::_communication.get() != NULL);
	  assertion(utils::MasterSlave::_communication->isConnected());

	  for(int i = 0; i < _oldInvJacobian.rows(); i++){
		  int rank = 0;
		  // find rank of processor that stores the result
		  // the second while is necessary if processors with no vertices are present
		  while(i >= _dimOffsets[rank+1]) rank++;
		  //if(rank > 0){ while(_dimOffsets[rank] == _dimOffsets[rank-1]) rank--;}

		  for(int j = 0; j < _qrV.cols(); j++){
			  // as we want to move to Eigen, copy
			  DataValues Jrow(_oldInvJacobian.cols(), 0.0);
			  for(int s = 0; s < _oldInvJacobian.cols(); s++){
				  if(!_hasNodesOnInterface) std::cout<<"should not be entered: i: "<<i<<", s: "<<s<<std::endl;
				  Jrow(s) = _oldInvJacobian(i,s);
			  }

			  // TODO: better: implement a reduce-operation (no loop over all slaves)
			  double res_ij = utils::MasterSlave::dot(Jrow, _matrixV.column(j));

			  // find proc that needs to store the result.
			  int local_row;
			  if(utils::MasterSlave::_rank == rank)
			  {
				  if(!_hasNodesOnInterface) std::cout<<"should not be entered, wrong rang found, rank: "<<rank<<std::endl;
				  local_row = i - _dimOffsets[rank];
				  tmpMatrix(local_row, j) = res_ij;
			  }
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
  preciceDebug("W_til * Z");
  assertion2(tmpMatrix.cols() == Z.rows(), tmpMatrix.cols(), Z.rows());
  if (not utils::MasterSlave::_masterMode && not utils::MasterSlave::_slaveMode) {
	  multiply(tmpMatrix, Z, _invJacobian);

  }else{

	  assertion(utils::MasterSlave::_communication.get() != NULL);
	  assertion(utils::MasterSlave::_communication->isConnected());

	  assertion(utils::MasterSlave::_cyclicCommLeft.get() != NULL);
	  assertion(utils::MasterSlave::_cyclicCommLeft->isConnected());

	  assertion(utils::MasterSlave::_cyclicCommRight.get() != NULL);
	  assertion(utils::MasterSlave::_cyclicCommRight->isConnected());


	  int nextProc = (utils::MasterSlave::_rank + 1) % utils::MasterSlave::_size;
	  int prevProc = (utils::MasterSlave::_rank -1 < 0) ? utils::MasterSlave::_size-1 : utils::MasterSlave::_rank -1;
	  int rows_rcv = (prevProc > 0) ? _dimOffsets[prevProc+1] - _dimOffsets[prevProc] : _dimOffsets[1];
	  //Matrix Wtil_rcv(rows_rcv, tmpMatrix.cols(),0.0);
	  Matrix Wtil_rcv(rows_rcv, getLSSystemCols(),0.0);

	  com::Request::SharedPointer requestSend;
	  com::Request::SharedPointer requestRcv;

	  // initiate asynchronous send operation of tmpMatrix --> nextProc (this data is needed in cycle 1)    dim: n_local x cols
	  if(tmpMatrix.size() > 0)
		  requestSend = utils::MasterSlave::_cyclicCommRight->aSend(&tmpMatrix(0,0), tmpMatrix.size(), 0);

	  // initiate asynchronous receive operation for W_til from previous processor --> W_til      dim: rows_rcv x cols
	  if(Wtil_rcv.size() > 0)
		  requestRcv = utils::MasterSlave::_cyclicCommLeft->aReceive(&Wtil_rcv(0,0), rows_rcv * tmpMatrix.cols(), 0);

	  // compute diagonal blocks where all data is local and no communication is needed
	  // compute block matrices of J_inv of size (n_til x n_til), n_til = local n
	  Matrix diagBlock(_matrixV.rows(),_matrixV.rows(), 0.0);
	  multiply(tmpMatrix, Z, diagBlock);
	  if(!_hasNodesOnInterface){
		  std::cout<<"diag block on empty proc in cycle 0. BLOCK: ("<<diagBlock.rows()<<","<<diagBlock.cols()<<") Z: ("<<Z.rows()<<","<<Z.cols()<<") Wtil: ("<<tmpMatrix.rows()<<","<<tmpMatrix.cols()<<")"<<std::endl;
	  }
	  // set block at corresponding row-index on proc
	  int off = _dimOffsets[utils::MasterSlave::_rank];
	  assertion2(_invJacobian.cols() == diagBlock.cols(), _invJacobian.cols(), diagBlock.cols());
	  for(int q = 0; q < diagBlock.rows(); q++)
		  for(int p = 0; p < _invJacobian.cols(); p++)
		  {
			  _invJacobian(q+off,p) = diagBlock(q,p);
		  }

	  /**
	   * cyclic send-receive operation
	   */
	  for(int cycle = 1; cycle < utils::MasterSlave::_size; cycle++){

		  // wait until W_til from previous processor is fully received
		  if(requestSend != NULL)
			  requestSend->wait();
		  if(requestRcv != NULL)
			  requestRcv->wait();


		  // Wtil_rcv is available - needed for local multiplication and hand over to next proc
		  Matrix Wtil_copy(Wtil_rcv);

		  // initiate async send to hand over Wtil to the next proc (this data will be needed in the next cycle)    dim: n_local x cols
		  if(cycle < utils::MasterSlave::_size-1){
			  if(Wtil_copy.size() > 0)
				  requestSend = utils::MasterSlave::_cyclicCommRight->aSend(&Wtil_copy(0,0), Wtil_copy.size(), 0);
			  //std::cout<<"proc["<<utils::MasterSlave::_rank<<"] sends to ["<<nextProc<<"] double* of size "<<Wtil_copy.size()<<" for the next cycle"<<std::endl;
		  }

		  // compute proc that owned Wtil_rcv at the very beginning for each cylce
		  int sourceProc_nextCycle = (utils::MasterSlave::_rank - (cycle+1) < 0) ?
				  utils::MasterSlave::_size + (utils::MasterSlave::_rank - (cycle+1)) : utils::MasterSlave::_rank - (cycle+1);
		  int sourceProc = (utils::MasterSlave::_rank - cycle < 0) ?
		  				  utils::MasterSlave::_size + (utils::MasterSlave::_rank - cycle) : utils::MasterSlave::_rank - cycle;
		  int rows_rcv_nextCycle = (sourceProc_nextCycle > 0) ? _dimOffsets[sourceProc_nextCycle+1] - _dimOffsets[sourceProc_nextCycle] : _dimOffsets[1];
		  rows_rcv = (sourceProc > 0) ? _dimOffsets[sourceProc+1] - _dimOffsets[sourceProc] : _dimOffsets[1];
		  Wtil_rcv = Matrix(rows_rcv_nextCycle, getLSSystemCols(),0.0);


		  // initiate asynchronous receive operation for W_til from previous processor --> W_til (this data is needed in the next cycle)
		  if(cycle < utils::MasterSlave::_size-1){
			  if(Wtil_rcv.size() > 0)
				  requestRcv = utils::MasterSlave::_cyclicCommLeft->aReceive(&Wtil_rcv(0,0), Wtil_rcv.size(), 0);
			//  std::cout<<"proc["<<utils::MasterSlave::_rank<<"] receives from ["<<prevProc<<"] double* of size "
			//		   <<(rows_rcv_nextCycle*tmpMatrix.cols())<<" from dest ["<<sourceProc_nextCycle<<"]"<<" for the next cycle"<<std::endl;
		  }


		  // compute block with new local data
		  Matrix block(rows_rcv, Z.cols(), 0.0);
		  multiply(Wtil_copy, Z, block);

		  if(!_hasNodesOnInterface){
		  		  std::cout<<"diag block on empty proc in cycle "<<cycle<<". BLOCK: ("<<block.rows()<<","<<block.cols()<<") Z: ("<<Z.rows()<<","<<Z.cols()<<") Wtil: ("<<Wtil_copy.rows()<<","<<Wtil_copy.cols()<<")"<<std::endl;
		  	  }

		  // set block at corresponding index in J_inv
		  // the row-offset of the current block is determined by the proc that sends the part of the W_til matrix
		  // note: the direction and ordering of the cyclic sending operation is chosen s.t. the computed block is
		  //       local on the current processor (in J_inv).
		  off = _dimOffsets[sourceProc];
		  assertion2(_invJacobian.cols() == block.cols(), _invJacobian.cols(), block.cols());
		  for(int q = 0; q < block.rows(); q++)
			  for(int p = 0; p < _invJacobian.cols(); p++)
			  {
				  _invJacobian(q+off,p) = block(q,p);
			  }

		//  std::cout<<"proc["<<utils::MasterSlave::_rank<<"] computed block ["<<sourceProc<<"]["<<utils::MasterSlave::_rank
		//		   <<"] inserted at off= "<<off<<" frobenius norm: "<<tarch::la::frobeniusNorm(block)<<std::endl;
	  }
  }

  _invJacobian = _invJacobian + _oldInvJacobian;

  DataValues negRes(_residuals);
  negRes *= -1.;

  /*
   * solve delta_x = - J_inv*residuals
   */
  if (not utils::MasterSlave::_masterMode && not utils::MasterSlave::_slaveMode) {
	  multiply(_invJacobian, negRes, xUpdate);

  }else{

  	  assertion(utils::MasterSlave::_communication.get() != NULL);
  	  assertion(utils::MasterSlave::_communication->isConnected());

  	  for(int i = 0; i < _invJacobian.rows(); i++){

		  int rank = 0;
		  // find rank of processor that stores the result
		  // the second while is necessary if processors with no vertices are present
		  while(i >= _dimOffsets[rank+1]) rank++;
		 // if(rank > 0){ while(_dimOffsets[rank] - _dimOffsets[rank-1] == 0) rank--;}

		  // as we want to move to Eigen, copy
		  DataValues Jrow(_invJacobian.cols(), 0.0);
		  for (int s = 0; s < _invJacobian.cols(); s++) {
			  if(!_hasNodesOnInterface) std::cout<<"should not be entered: i: "<<i<<", s: "<<s<<std::endl;
			  Jrow(s) = _invJacobian(i,s);
		  }


		  // TODO: better: implement a reduce-operation (no loop over all slaves)
		  double up_ij = utils::MasterSlave::dot(Jrow, negRes);

		  // find proc that needs to store the result.
		  int local_row;
		  if(utils::MasterSlave::_rank == rank)
		  {
			  if(!_hasNodesOnInterface) std::cout<<"should not be entered, wrong rang found, rank: "<<rank<<std::endl;
			  local_row = i - _dimOffsets[rank];
			  xUpdate(local_row) = up_ij;
		  }
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
  // store inverse Jacobian
//  _matrixWriter.write(_invJacobian);
  _oldInvJacobian = _invJacobian;
}

}}} // namespace precice, cplscheme, impl
