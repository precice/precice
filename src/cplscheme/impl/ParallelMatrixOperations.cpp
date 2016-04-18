/*
 * ParallelMatrixOperations.cpp
 *
 *  Created on: Aug 21, 2015
 *      Author: Klaudius Scheufele
 */
// Copyright (C) 2015 Universität Stuttgart
#ifndef PRECICE_NO_MPI

#include "ParallelMatrixOperations.hpp"
#include "utils/Globals.hpp"
#include "tarch/la/MatrixVectorOperations.h"
#include "tarch/la/MatrixMatrixOperations.h"
#include "tarch/la/TransposedMatrix.h"
#include "utils/Dimensions.hpp"
#include "tarch/la/Scalar.h"
#include "utils/MasterSlave.hpp"
#include "com/MPIPortsCommunication.hpp"
#include "com/Communication.hpp"
#include "logging/Logger.hpp"
#include "Eigen/Dense"


namespace precice {
namespace cplscheme {
namespace impl {

logging::Logger ParallelMatrixOperations::
      _log("precice::cplscheme::impl::ParallelMatrixOperations");


ParallelMatrixOperations::ParallelMatrixOperations() :
_cyclicCommLeft(nullptr),
_cyclicCommRight(nullptr),
_needCycliclComm(true)
{}

void ParallelMatrixOperations::initialize(
		com::Communication::SharedPointer leftComm,
		com::Communication::SharedPointer rightComm,
		bool needCyclicComm)
{
	preciceTrace("initialize()");

	_needCycliclComm = needCyclicComm;
	if(utils::MasterSlave::_masterMode ||utils::MasterSlave::_slaveMode){

		_cyclicCommLeft = leftComm;
		_cyclicCommRight = rightComm;

		if(_needCycliclComm){
		  assertion(_cyclicCommLeft.get() != NULL); assertion(_cyclicCommLeft->isConnected());
		  assertion(_cyclicCommRight.get() != NULL); assertion(_cyclicCommRight->isConnected());
		}
	}
}




// =============================== TARCH METHODS (NOT NEEDED ANY  MORE) =======================

void ParallelMatrixOperations::multiply
(
	TarchMatrix& leftMatrix,
	TarchMatrix& rightMatrix,
	TarchMatrix& result,
	const std::vector<int>& offsets,
	int p, int q, int r)
{
	preciceTrace("multiply()");
	assertion(result.cols() == rightMatrix.cols(), result.cols(), rightMatrix.cols());
	assertion(leftMatrix.cols() == rightMatrix.rows(), leftMatrix.cols(), rightMatrix.rows());

	// if serial computation on single processor, i.e, no master-slave mode
	if( not utils::MasterSlave::_masterMode && not utils::MasterSlave::_slaveMode){
		tarch::la::multiply(leftMatrix, rightMatrix, result);

	// if parallel computation on p processors, i.e., master-slave mode
	}else{
		assertion(utils::MasterSlave::_communication.get() != NULL);
		assertion(utils::MasterSlave::_communication->isConnected());

		// The result matrix is of size (p x r)
		// if p equals r (and p = global_n), we have to perform the
		// cyclic communication with block-wise matrix-matrix multiplication
		if(p == r){
		  assertion(_needCycliclComm);
			assertion(_cyclicCommLeft.get() != NULL); assertion(_cyclicCommLeft->isConnected());
			assertion(_cyclicCommRight.get() != NULL); assertion(_cyclicCommRight->isConnected());

			_multiplyNN(leftMatrix, rightMatrix, result, offsets, p, q, r);

		// case p != r, i.e., usually p = number of columns of the least squares system
		// perform parallel multiplication based on dot-product
		}else{
			_multiplyNM(leftMatrix, rightMatrix, result, offsets, p, q, r);
		}
	}
}

void ParallelMatrixOperations::multiply
(
	TarchMatrix& leftMatrix,
	TarchColumnMatrix& rightMatrix,
	TarchMatrix& result,
	const std::vector<int>& offsets,
	int p, int q, int r)
{
	preciceTrace("multiply()");
	TarchMatrix rM(rightMatrix.rows(), rightMatrix.cols());
	for(int i = 0; i < rightMatrix.rows(); i++)
		for(int j = 0; j < rightMatrix.cols(); j++){
			rM(i,j) = rightMatrix(i,j);
		}
	multiply(leftMatrix, rM, result, offsets, p, q, r);
}

// @brief multiplies matrices based on a dot-product computation with a rectangular result matrix
void ParallelMatrixOperations::_multiplyNM(
		TarchMatrix& leftMatrix, TarchMatrix& rightMatrix, TarchMatrix& result, const std::vector<int>& offsets, int p, int q, int r)
{
	preciceTrace("multiplyNM()");
	for(int i = 0; i < leftMatrix.rows(); i++){
	  int rank = 0;
	  // find rank of processor that stores the result
	  // Note: the >'=' here is crucial: In case some procs do not have any vertices,
	  // this while loop continues incrementing rank if entries in offsets are equal, i.e.,
	  // it runs to the next non-empty proc.
	  while(i >= offsets[rank+1]) rank++;

	  TarchVector lMRow(leftMatrix.cols(), 0.0);
	  for(int s = 0; s < leftMatrix.cols(); s++){
		  lMRow(s) = leftMatrix(i,s);
	  }

	  for(int j = 0; j < r; j++){

		  TarchVector rMCol(rightMatrix.rows(), 0.0);
		  for(int s = 0; s < rightMatrix.rows(); s++){
			  rMCol(s) = rightMatrix(s,j);
		  }
		  double res_ij = utils::MasterSlave::dot(lMRow, rMCol);

		  // find proc that needs to store the result.
		  int local_row;
		  if(utils::MasterSlave::_rank == rank)
		  {
			  local_row = i - offsets[rank];
			  result(local_row, j) = res_ij;
		  }
	  }
	}
}


// @brief multiplies matrices based on a cyclic communication and block-wise matrix multiplication with a quadratic result matrix
void ParallelMatrixOperations::_multiplyNN(
		TarchMatrix& leftMatrix, TarchMatrix& rightMatrix, TarchMatrix& result, const std::vector<int>& offsets, int p, int q, int r)
{
	preciceTrace("multiplyNN()");

	assertion(_needCycliclComm);
	/*
	 * For multiplication W_til * Z = J
	 * -----------------------------------------------------------------------
	 * p = r = n_global, q = m
	 *
	 * leftMatrix:  local: (n_local x m) 		global: (n_global x m)
	 * rightMatrix: local: (m x n_local) 		global: (m x n_global)
	 * result: 		local: (n_global x n_local) global: (n_global x n_global)
	 * -----------------------------------------------------------------------
	 */


	//int nextProc = (utils::MasterSlave::_rank + 1) % utils::MasterSlave::_size;
	int prevProc = (utils::MasterSlave::_rank -1 < 0) ? utils::MasterSlave::_size-1 : utils::MasterSlave::_rank -1;
	int rows_rcv = (prevProc > 0) ? offsets[prevProc+1] - offsets[prevProc] : offsets[1];
	TarchMatrix leftMatrix_rcv(rows_rcv, q, 0.0);

	com::Request::SharedPointer requestSend;
	com::Request::SharedPointer requestRcv;

	// initiate asynchronous send operation of leftMatrix (W_til) --> nextProc (this data is needed in cycle 1)    dim: n_local x cols
	if(leftMatrix.size() > 0)
		requestSend = _cyclicCommRight->aSend(&leftMatrix(0,0), leftMatrix.size(), 0);

	// initiate asynchronous receive operation for leftMatrix (W_til) from previous processor --> W_til      dim: rows_rcv x cols
	if(leftMatrix_rcv.size() > 0)
		requestRcv = _cyclicCommLeft->aReceive(&leftMatrix_rcv(0,0), rows_rcv * leftMatrix.cols(), 0);

	// compute diagonal blocks where all data is local and no communication is needed
	// compute block matrices of J_inv of size (n_til x n_til), n_til = local n
	TarchMatrix diagBlock(leftMatrix.rows(), leftMatrix.rows(), 0.0);
	tarch::la::multiply(leftMatrix, rightMatrix, diagBlock);

	// set block at corresponding row-index on proc
	int off = offsets[utils::MasterSlave::_rank];
	assertion(result.cols() == diagBlock.cols(), result.cols(), diagBlock.cols());
	for(int ii = 0; ii < diagBlock.rows(); ii++)
		for(int jj = 0; jj < result.cols(); jj++)
		{
		  result(ii+off, jj) = diagBlock(ii, jj);
		}

	/**
	 * cyclic send-receive operation
	 */
	for(int cycle = 1; cycle < utils::MasterSlave::_size; cycle++){

		// wait until W_til from previous processor is fully received
		if(requestSend != NULL) requestSend->wait();
		if(requestRcv != NULL)  requestRcv->wait();
	
		// leftMatrix (leftMatrix_rcv) is available - needed for local multiplication and hand over to next proc
		TarchMatrix leftMatrix_copy(leftMatrix_rcv);
	
		// initiate async send to hand over leftMatrix (W_til) to the next proc (this data will be needed in the next cycle)    dim: n_local x cols
		if(cycle < utils::MasterSlave::_size-1){
		  if(leftMatrix_copy.size() > 0)
			  requestSend = _cyclicCommRight->aSend(&leftMatrix_copy(0,0), leftMatrix_copy.size(), 0);
		}
	
		// compute proc that owned leftMatrix_rcv (Wtil_rcv) at the very beginning for each cylce
		int sourceProc_nextCycle = (utils::MasterSlave::_rank - (cycle+1) < 0) ?
			  utils::MasterSlave::_size + (utils::MasterSlave::_rank - (cycle+1)) : utils::MasterSlave::_rank - (cycle+1);
	
		int sourceProc = (utils::MasterSlave::_rank - cycle < 0) ?
			  utils::MasterSlave::_size + (utils::MasterSlave::_rank - cycle) : utils::MasterSlave::_rank - cycle;
	
		int rows_rcv_nextCycle = (sourceProc_nextCycle > 0) ? offsets[sourceProc_nextCycle+1] - offsets[sourceProc_nextCycle] : offsets[1];
		rows_rcv = (sourceProc > 0) ? offsets[sourceProc+1] - offsets[sourceProc] : offsets[1];
		leftMatrix_rcv = TarchMatrix(rows_rcv_nextCycle, q, 0.0);
	
	
		// initiate asynchronous receive operation for leftMatrix (W_til) from previous processor --> W_til (this data is needed in the next cycle)
		if(cycle < utils::MasterSlave::_size-1){
		  if(leftMatrix_rcv.size() > 0) // only receive data, if data has been sent
			  requestRcv = _cyclicCommLeft->aReceive(&leftMatrix_rcv(0,0), leftMatrix_rcv.size(), 0);
		}
	
		// compute block with new local data
		TarchMatrix block(rows_rcv, rightMatrix.cols(), 0.0);
		tarch::la::multiply(leftMatrix_copy, rightMatrix, block);
	
		// set block at corresponding index in J_inv
		// the row-offset of the current block is determined by the proc that sends the part of the W_til matrix
		// note: the direction and ordering of the cyclic sending operation is chosen s.t. the computed block is
		//       local on the current processor (in J_inv).
		off = offsets[sourceProc];
		assertion(result.cols() == block.cols(), result.cols(), block.cols());
		for(int ii = 0; ii < block.rows(); ii++)
		  for(int jj = 0; jj < result.cols(); jj++)
		  {
			  result(ii+off, jj) = block(ii, jj);
		  }
	}
}


void ParallelMatrixOperations::multiply
(
	TarchMatrix& A,
	TarchVector& v,
	TarchVector& result,
	const std::vector<int>& offsets,
	int p, int q)
{
	preciceTrace("multiply()");
	assertion(v.size() == A.cols(), v.size(), A.cols());

	// if serial computation on single processor, i.e, no master-slave mode
	if( not utils::MasterSlave::_masterMode && not utils::MasterSlave::_slaveMode){
		tarch::la::multiply(A, v, result);

	// if parallel computation on p processors, i.e., master-slave mode
	}else{
		assertion(utils::MasterSlave::_communication.get() != NULL);
		assertion(utils::MasterSlave::_communication->isConnected());


	  	  for(int i = 0; i < A.rows(); i++){
			  int rank = 0;
			  // find rank of processor that stores the result
			  // the second while is necessary if processors with no vertices are present
			  while(i >= offsets[rank+1]) rank++;

			  TarchVector Arow(A.cols(), 0.0);
			  for (int s = 0; s < A.cols(); s++) {
				  Arow(s) = A(i,s);
			  }
			  double up_ij = utils::MasterSlave::dot(Arow, v);

			  // find proc that needs to store the result.
			  int local_row;
			  if(utils::MasterSlave::_rank == rank)
			  {
				  local_row = i - offsets[rank];
				  result(local_row) = up_ij;
			  }
	  	  }
	}
}


}}} // namespace precice, cplscheme, impl

#endif
