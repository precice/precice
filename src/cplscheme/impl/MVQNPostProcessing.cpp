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
  _invJ(),
  _oldInvJ(),
  _V(),
  _W(),
  _Wtil(),
  _Z(),
  _dimOffsets()
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

		assertion(utils::MasterSlave::_communication.get() != NULL);
		assertion(utils::MasterSlave::_communication->isConnected());


		/**
		 *  make dimensions public to all procs,
		 *  last entry _dimOffsets[MasterSlave::_size] holds the global dimension, global,n
		 */
		_dimOffsets.resize(utils::MasterSlave::_size + 1);
		if (utils::MasterSlave::_slaveMode) {
			utils::MasterSlave::_communication->send(entries, 0);
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
		}
		global_n = _dimOffsets.back();
	}
  
  _invJacobian = Matrix(global_n, entries, init);
  _oldInvJacobian = Matrix(global_n, entries, init);

  // debug
  //if(utils::MasterSlave::_masterMode){
	  _invJ = Matrix(global_n, global_n, init);
	  _oldInvJ = Matrix(global_n, global_n, init);
  //}
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

  //DataMatrix Z;
  Matrix Z(_matrixV.cols(), _matrixV.rows(), 0.0);
  bool linearDependence = true;
  
  while (linearDependence) {
		linearDependence = false;
		//Z.clear();

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
				for(int p = 0; p < __ytmpVec.size(); p++)
					Z(p,i) = __ytmpVec(p);
				//Z.append(__ytmpVec);
				__matrixQRow.clear();
			}
		}
	}


  //debug
  DataValues _R(_invJacobian.rows(), 0.0);


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

	  // debug ################################################################################### debug
	  // #######
	  if(utils::MasterSlave::_slaveMode){
		  Matrix send_matrixV(_matrixV.rows(), _matrixV.cols(), 0.0);
		  for (int i = 0; i < _matrixV.rows(); i++)
			  for (int j = 0; j < _matrixV.cols(); j++) {
				  send_matrixV(i, j) = _matrixV(i, j);
			  }
		  Matrix send_matrixW(_matrixW.rows(), _matrixW.cols(), 0.0);
		  for (int i = 0; i < _matrixW.rows(); i++)
			  for (int j = 0; j < _matrixW.cols(); j++) {
				  send_matrixW(i, j) = _matrixW(i, j);
			  }
//		  utils::MasterSlave::_communication->send(&_invJacobian(0,0), _invJacobian.size(), 0);
		  utils::MasterSlave::_communication->send(&send_matrixV(0,0), send_matrixV.size(), 0);
		  utils::MasterSlave::_communication->send(&send_matrixW(0,0), send_matrixW.size(), 0);
		  utils::MasterSlave::_communication->send(&Z(0,0), Z.size(), 0);
		  utils::MasterSlave::_communication->send(&_residuals(0), _residuals.size(), 0);
	  }

	  _V = Matrix(_invJacobian.rows(), _matrixV.cols(), 0.0);
	  _W = Matrix(_invJacobian.rows(), _matrixW.cols(), 0.0);
	  _Z = Matrix(_matrixV.cols(), _invJacobian.rows(), 0.0);
	  _Wtil = Matrix(_invJacobian.rows(), _matrixW.cols(), 0.0);

	  if(utils::MasterSlave::_masterMode){

		  for(int p = 0; p < _matrixV.rows(); p++)
			  for(int q = 0; q < _matrixV.cols(); q++)
			  {
				  _V(p,q) = _matrixV(p,q);
			  }
		  for(int p = 0; p < _matrixW.rows(); p++)
			  for(int q = 0; q < _matrixW.cols(); q++)
			  {
				  _W(p,q) = _matrixW(p,q);
			  }
		  for(int p = 0; p < Z.rows(); p++)
			  for(int q = 0; q < Z.cols(); q++)
			  {
				  _Z(p,q) = Z(p,q);
			  }
//		  for(int p = 0; p < _invJacobian.rows(); p++)
//			  for(int q = 0; q < _invJacobian.cols(); q++)
//			  {
//				  _invJ(p,q) = _invJacobian(p,q);
//			  }
		  for(int p = 0; p < _residuals.size(); p++)
			  _R(p) = _residuals(p);

		  for(int rankSlave = 1; rankSlave < utils::MasterSlave::_size; rankSlave++){
			  int r = _dimOffsets[rankSlave+1]-_dimOffsets[rankSlave];
//			  Matrix J(_invJacobian.rows(), r);
//			  std::cout<<"receive J from ["<<rankSlave<<"], size="<<J.size()<<", r="<<r<<std::endl;
//			  utils::MasterSlave::_communication->receive(&J(0,0), J.size(), rankSlave);
			  int o = _dimOffsets[rankSlave];
//		      for(int p = 0; p < J.rows(); p++)
//		    	  for(int q = 0; q < J.cols(); q++)
//		    	  {
//		    		  _invJ(p,q+o) = J(p,q);
//		    	  }

		      Matrix V(r, _matrixV.cols());
		      std::cout<<"receive V from ["<<rankSlave<<"], size="<<V.size()<<", r="<<r<<std::endl;
		      utils::MasterSlave::_communication->receive(&V(0,0), V.size(), rankSlave);
		      for(int p = 0; p < V.rows(); p++)
				  for(int q = 0; q < V.cols(); q++)
				  {
					  _V(p+o,q) = V(p,q);
				  }

		      Matrix W(r, _matrixW.cols());
		      std::cout<<"receive W from ["<<rankSlave<<"], size="<<W.size()<<", r="<<r<<std::endl;
			  utils::MasterSlave::_communication->receive(&W(0,0), W.size(), rankSlave);
			  for(int p = 0; p < W.rows(); p++)
				  for(int q = 0; q < W.cols(); q++)
				  {
					  _W(p+o,q) = W(p,q);
				  }

			  Matrix zz(_matrixV.cols(), r);
			  std::cout<<"receive Z from ["<<rankSlave<<"], size="<<zz.size()<<", r="<<r<<std::endl;
			  utils::MasterSlave::_communication->receive(&zz(0,0), zz.size(), rankSlave);
			  for(int p = 0; p < zz.rows(); p++)
				  for(int q = 0; q < zz.cols(); q++)
				  {
					  _Z(p,q+o) = zz(p,q);
				  }

			  DataValues res(r, 0.0);
			  std::cout<<"receive res from ["<<rankSlave<<"], size="<<res.size()<<", r="<<r<<std::endl;
			  utils::MasterSlave::_communication->receive(&res(0), res.size(), rankSlave);
			  for(int p = 0; p < res.size(); p++)
				  _R(p+o) = res(p);
		  }

		  multiply(_oldInvJ, _V, _Wtil);

		  std::cout<<"  _Z from server: "<<tarch::la::frobeniusNorm(_Z)<<std::endl;
	  }
	  utils::MasterSlave::broadcast(&_Wtil(0,0), _Wtil.size());
	  // #######
	  // debug ################################################################################### debug



	  for(int i = 0; i < _oldInvJacobian.rows(); i++){
		  // find rank of processor that stores the result
		  int rank = 0;
		  while(i >= _dimOffsets[rank+1]) rank++;
		  for(int j = 0; j < _matrixV.cols(); j++){
			  // as we want to move to Eigen, copy
			  DataValues Jrow(_oldInvJacobian.cols(), 0.0);
			  for(int s = 0; s < _oldInvJacobian.cols(); s++)
				  Jrow(s) = _oldInvJacobian(i,s);

			  // TODO: better: implement a reduce-operation (no loop over all slaves)
			  double res_ij = utils::MasterSlave::dot(Jrow, _matrixV.column(j));

			  // find proc that needs to store the result.
			  int local_row;
			  if(utils::MasterSlave::_rank == rank)
			  {
				  local_row = i - _dimOffsets[rank];
				  tmpMatrix(local_row, j) = res_ij;
			  }
		  }
	  }

  }

  // tmpMatrix = (W-J_inv_n*V)
  tmpMatrix *= -1.;
  tmpMatrix = tmpMatrix + _matrixW;



  // debug ################################################################################### debug
  	  // #######
  	  // VALIDATE

	  if(utils::MasterSlave::_masterMode)
	  {
		  _Wtil *= -1;
		  _Wtil = _Wtil + _W;
	  }
	  utils::MasterSlave::broadcast(&_Wtil(0,0), _Wtil.size());
	  std::cout<<"proc["<<utils::MasterSlave::_rank<<"]  _Wtil from server: "<<tarch::la::frobeniusNorm(_Wtil)<<std::endl;

  	  bool failed = false;
        double largestDiff = 0.;
        double diff = 0., val1 = 0., val2 = 0.;
  	  int o = _dimOffsets[utils::MasterSlave::_rank];
  	  for(int p = 0; p < tmpMatrix.rows(); p++)
  		  for(int q = 0; q < tmpMatrix.cols(); q++)
  		  {
  				if(!tarch::la::equals(_Wtil(p+o,q), tmpMatrix(p,q), 1e-12)) {
  					failed = true;
  					diff = _Wtil(p,q+o) - tmpMatrix(p,q);
  					if (largestDiff < diff) {
  						largestDiff = diff;
  						val1 = tmpMatrix(p,q);
  						val2 = _Wtil(p,q+o);
  					}
  				}
  		  }
  	  if(failed)
  	  {
  		  std::cerr<<"proc["<<utils::MasterSlave::_rank<<"] validation failed for Wtil.\n  - master-slave: "<<val1<<"\n  - server: "<<val2<<"\n  - (max-)diff: "<<largestDiff<<std::endl;
  		  _infostream<<"proc["<<utils::MasterSlave::_rank<<"] validation failed for Wtil.\n  - master-slave: "<<val1<<"\n  - server: "<<val2<<"\n  - (max-)diff: "<<largestDiff<<std::endl;
  	  }
  	  // #######
  	  // debug ################################################################################### debug
  
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

	  assertion(utils::MasterSlave::_cyclicCommLeft.get() != NULL);
	  assertion(utils::MasterSlave::_cyclicCommLeft->isConnected());

	  assertion(utils::MasterSlave::_cyclicCommRight.get() != NULL);
	  assertion(utils::MasterSlave::_cyclicCommRight->isConnected());


	  int nextProc = (utils::MasterSlave::_rank + 1) % utils::MasterSlave::_size;
	  int prevProc = (utils::MasterSlave::_rank -1 < 0) ? utils::MasterSlave::_size-1 : utils::MasterSlave::_rank -1;
	  int rows_rcv = (prevProc > 0) ? _dimOffsets[prevProc+1] - _dimOffsets[prevProc] : _dimOffsets[1];
	  Matrix Wtil_rcv(rows_rcv, tmpMatrix.cols(),0.0);

	  // initiate asynchronous send operation of tmpMatrix --> nextProc (this data is needed in cycle 1)    dim: n_local x cols
	  com::Request::SharedPointer requestSend = utils::MasterSlave::_cyclicCommRight->aSend(&tmpMatrix(0,0), tmpMatrix.size(), 0);

	  // initiate asynchronous receive operation for W_til from previous processor --> W_til      dim: rows_rcv x cols
	  com::Request::SharedPointer requestRcv = utils::MasterSlave::_cyclicCommLeft->aReceive(&Wtil_rcv(0,0), rows_rcv * tmpMatrix.cols(), 0);

	  // compute diagonal blocks where all data is local and no communication is needed
	  // compute block matrices of J_inv of size (n_til x n_til), n_til = local n
	  Matrix diagBlock(_matrixV.rows(),_matrixV.rows(), 0.0);
	  multiply(tmpMatrix, Z, diagBlock);
	  // set block at corresponding row-index on proc
	  int off = _dimOffsets[utils::MasterSlave::_rank];
	  assertion2(_invJacobian.cols() == diagBlock.cols(), _invJacobian.cols(), diagBlock.cols());
	  for(int q = 0; q < diagBlock.rows(); q++)
		  for(int p = 0; p < _invJacobian.cols(); p++)
		  {
			  _invJacobian(q+off,p) = diagBlock(q,p);
		  }

	  std::cout<<"computed block ["<<utils::MasterSlave::_rank<<"]["<<utils::MasterSlave::_rank<<"] inserted at off= "<<off<<" frobenius norm: "<<tarch::la::frobeniusNorm(diagBlock)<<std::endl;


	  // cyclic send-receive operation
	  for(int cycle = 1; cycle < utils::MasterSlave::_size; cycle++){

		  // wait until W_til from previous processor is fully received
		  requestSend->wait();
		  requestRcv->wait();

		  //debug
		  int sourceProc = (utils::MasterSlave::_rank - cycle < 0) ?
		  		  				  utils::MasterSlave::_size + (utils::MasterSlave::_rank - cycle) : utils::MasterSlave::_rank - cycle;
		  std::cout<<"proc["<<utils::MasterSlave::_rank<<"] receives Wtil from ["<<sourceProc<<"]: "<<tarch::la::frobeniusNorm(Wtil_rcv)<<" size: "<<Wtil_rcv.rows()<<std::endl;

		  // Wtil_rcv is available - needed for local multiplication and hand over to next proc
		  Matrix Wtil_copy(Wtil_rcv.rows(), Wtil_rcv.cols(), 0.0);
		  for(int ii = 0; ii < Wtil_rcv.rows(); ii++)
			  for(int jj = 0; jj < Wtil_rcv.cols(); jj++)
				  Wtil_copy(ii,jj) = Wtil_rcv(ii,jj);


/*
		  {
				      int i = 0;
				      char hostname[256];
				      gethostname(hostname, sizeof(hostname));
				      printf("PID %d on %s ready for attach\n", getpid(), hostname);
				      fflush(stdout);
				      while (0 == i)
				          sleep(5);
		  }
*/
		  // initiate async send to hand over Wtil to the next proc (this data will be needed in the next cycle)    dim: n_local x cols
		  if(cycle < utils::MasterSlave::_size-1){
			  requestSend = utils::MasterSlave::_cyclicCommRight->aSend(&Wtil_copy(0,0), Wtil_copy.size(), 0);
			  std::cout<<"proc["<<utils::MasterSlave::_rank<<"] sends to ["<<nextProc<<"] double* of size "<<Wtil_copy.size()<<" for the next cycle"<<std::endl;
		  }

		  // compute proc that owned Wtil_rcv at the very beginning for each cylce
		  int sourceProc_nextCycle = (utils::MasterSlave::_rank - (cycle+1) < 0) ?
				  utils::MasterSlave::_size + (utils::MasterSlave::_rank - (cycle+1)) : utils::MasterSlave::_rank - (cycle+1);
		  sourceProc = (utils::MasterSlave::_rank - cycle < 0) ?
		  				  utils::MasterSlave::_size + (utils::MasterSlave::_rank - cycle) : utils::MasterSlave::_rank - cycle;
		  int rows_rcv_nextCycle = (sourceProc_nextCycle > 0) ? _dimOffsets[sourceProc_nextCycle+1] - _dimOffsets[sourceProc_nextCycle] : _dimOffsets[1];
		  rows_rcv = (sourceProc > 0) ? _dimOffsets[sourceProc+1] - _dimOffsets[sourceProc] : _dimOffsets[1];
		  Wtil_rcv = Matrix(rows_rcv_nextCycle, tmpMatrix.cols(),0.0);


		  // initiate asynchronous receive operation for W_til from previous processor --> W_til (this data is needed in the next cycle)
		  if(cycle < utils::MasterSlave::_size-1){
			  requestRcv = utils::MasterSlave::_cyclicCommLeft->aReceive(&Wtil_rcv(0,0), rows_rcv_nextCycle * tmpMatrix.cols(), 0);
			  std::cout<<"proc["<<utils::MasterSlave::_rank<<"] receives from ["<<prevProc<<"] double* of size "
					   <<(rows_rcv_nextCycle*tmpMatrix.cols())<<" from dest ["<<sourceProc_nextCycle<<"]"<<" for the next cycle"<<std::endl;
		  }


		  // compute block with new local data
		  Matrix block(rows_rcv, Z.cols(), 0.0);
		  multiply(Wtil_copy, Z, block);

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

		  std::cout<<"proc["<<utils::MasterSlave::_rank<<"] computed block ["<<sourceProc<<"]["<<utils::MasterSlave::_rank
				   <<"] inserted at off= "<<off<<" frobenius norm: "<<tarch::la::frobeniusNorm(block)<<std::endl;
	  }

	  // debug ################################################################################### debug
	  // #######


	  if(utils::MasterSlave::_masterMode)
	  {
		  multiply(_Wtil, _Z, _invJ);
		  std::cout<<"  _invJ: "<<tarch::la::frobeniusNorm(_invJ)<<std::endl;
		  std::cout<<"  _Z: "<<tarch::la::frobeniusNorm(_Z)<<std::endl;
		  std::cout<<"  _Wtil: "<<tarch::la::frobeniusNorm(_Wtil)<<std::endl;
	  }
	  utils::MasterSlave::broadcast(&_invJ(0,0), _invJ.size());


	  std::cout<<"proc["<<utils::MasterSlave::_rank<<"] _invJacobian: "<<tarch::la::frobeniusNorm(_invJacobian)<<std::endl;
	  Matrix partJacobian(_invJacobian.rows(), _invJacobian.cols());
	  for(int q = 0; q < partJacobian.rows(); q++)
		  for(int p = 0; p < partJacobian.cols(); p++)
		  {
			  partJacobian(q,p) = _invJ(q,p+(_dimOffsets[utils::MasterSlave::_rank]));
		  }

	  std::cout<<"proc["<<utils::MasterSlave::_rank<<"] _invJacobian from server: "<<tarch::la::frobeniusNorm(partJacobian)<<std::endl;

	  for(int i = 0; i < utils::MasterSlave::_size; i++){
		  int oo = _dimOffsets[i];
		  Matrix bb(_dimOffsets[i+1]-oo, partJacobian.cols(), 0.0);
		  for(int q = 0; q < bb.rows(); q++)
			  for(int p = 0; p < partJacobian.cols(); p++)
			  {
				  bb(q,p) = partJacobian(q+oo,p);
			  }
		  std::cout<<"proc["<<utils::MasterSlave::_rank<<"] block ["<<i<<"]["<<utils::MasterSlave::_rank<<"] from server: "<<tarch::la::frobeniusNorm(bb)<<std::endl;
	  }
	  	  // VALIDATE

	  	  bool failed = false;
	        double largestDiff = 0.;
	        double diff = 0., val1 = 0., val2 = 0.;
	  	  int o = _dimOffsets[utils::MasterSlave::_rank];
	  	  int ii = 0, jj= 0, count = 0;
	  	  for(int p = 0; p < _invJacobian.rows(); p++)
	  		  for(int q = 0; q < _invJacobian.cols(); q++)
	  		  {
	  				if(!tarch::la::equals(_invJ(p,q+o), _invJacobian(p,q), 1e-10)) {
	  					failed = true;
	  					count++;
	  					diff = _invJ(p,q+o) - _invJacobian(p,q);
	  					if (largestDiff < diff) {
	  						largestDiff = diff;
	  						val1 = _invJacobian(p,q);
	  						val2 = _invJ(p,q+o);
	  						ii = p; jj = q;
	  					}
	  				}
	  		  }
	  	  if(failed)
	  	  {
	  		  std::cerr<<"proc["<<utils::MasterSlave::_rank<<"] validation failed for _invJacobian at ("<<ii<<", "<<jj<<") or ("<<ii<<", "<<jj+o<<"), resp. Total: "<<count<<".\n  - master-slave: "<<val1<<"\n  - server: "<<val2<<"\n  - (max-)diff: "<<largestDiff<<std::endl;
	  		  _infostream<<"proc["<<utils::MasterSlave::_rank<<"] validation failed for _invJacobian at ("<<ii<<", "<<jj<<") or ("<<ii<<", "<<jj+o<<"), resp. Total: "<<count<<".\n  - master-slave: "<<val1<<"\n  - server: "<<val2<<"\n  - (max-)diff: "<<largestDiff<<std::endl;
	  	  }
	  	  if(utils::MasterSlave::_masterMode)
		  {
			  _invJ = _invJ + _oldInvJ;
		  }
	  	  utils::MasterSlave::broadcast(&_invJ(0,0), _invJ.size());
	  	// #######
	  	// debug ################################################################################### debug^

  }

  _invJacobian = _invJacobian + _oldInvJacobian;

  DataValues negRes(_residuals);
  negRes *= -1.;


  //debug
  DataValues _xUp(_R.size(), 0.0);


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

  	  std::cout<<"multiply J_inv * Res"<<std::endl;

  	  for(int i = 0; i < _invJacobian.rows(); i++){
  		  // find rank of processor that stores the result
		  int rank = 0;
		  while(i >= _dimOffsets[rank+1]) rank++;
		  // as we want to move to Eigen, copy
		  DataValues Jrow(_invJacobian.cols(), 0.0);
		  for (int s = 0; s < _invJacobian.cols(); s++) {
			  Jrow(s) = _invJacobian(i,s);
		  }
		  // TODO: better: implement a reduce-operation (no loop over all slaves)
		  double up_ij = utils::MasterSlave::dot(Jrow, negRes);

		  // find proc that needs to store the result.
		  int local_row;
		  if(utils::MasterSlave::_rank == rank)
		  {
			  local_row = i - _dimOffsets[rank];
			  xUpdate(local_row) = up_ij;
		  }
  	  }

  	  // debug ################################################################################### debug
	  // #######
	  if(utils::MasterSlave::_masterMode)
	  {
		  _R *= -1;
		  multiply(_invJ, _R, _xUp);
	  }
	  utils::MasterSlave::broadcast(&_xUp(0), _xUp.size());

	  // VALIDATE

	  	bool failed = false;
		double largestDiff = 0.;
		double diff = 0., val1 = 0., val2 = 0.;
	  int o = _dimOffsets[utils::MasterSlave::_rank];
	  for(int p = 0; p < xUpdate.size(); p++)
	  {
			if(!tarch::la::equals(_xUp(p+o), xUpdate(p), 1e-10)) {
				failed = true;
				diff = _xUp(p+o) - xUpdate(p);
				if (largestDiff < diff) {
					largestDiff = diff;
					val1 = xUpdate(p);
					val2 = _xUp(p+o);
				}
			}
	  }
	  if(failed)
	  {
		  std::cerr<<"proc["<<utils::MasterSlave::_rank<<"] validation failed for xUpdate.\n  - master-slave: "<<val1<<"\n  - server: "<<val2<<"\n  - (max-)diff: "<<largestDiff<<std::endl;
		  _infostream<<"proc["<<utils::MasterSlave::_rank<<"] validation failed for xUpdate.\n  - master-slave: "<<val1<<"\n  - server: "<<val2<<"\n  - (max-)diff: "<<largestDiff<<std::endl;
	  }


	  for(int i = 0; i < xUpdate.size(); i++){
		  xUpdate(i) = _xUp(i+o);
	  }


	  // #######
	  // debug ################################################################################### debug^

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
