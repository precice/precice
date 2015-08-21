/*
 * ParallelMatrixOperations.cpp
 *
 *  Created on: Aug 21, 2015
 *      Author: scheufks
 */
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
#include "tarch/logging/Log.h"
#include "Eigen/Dense"


namespace precice {
namespace cplscheme {
namespace impl {

tarch::logging::Log ParallelMatrixOperations::
      _log("precice::cplscheme::impl::ParallelMatrixOperations");


ParallelMatrixOperations::ParallelMatrixOperations() :
_cyclicCommLeft(nullptr),
_cyclicCommRight(nullptr)
{}

ParallelMatrixOperations::initialize(
		com::Communication::SharedPointer leftComm,
	   	com::Communication::SharedPointer rightComm)
{
	if(utils::MasterSlave::_masterMode ||utils::MasterSlave::_slaveMode){

		assertion(_cyclicCommLeft.get() != NULL); assertion(_cyclicCommLeft->isConnected());
		assertion(_cyclicCommRight.get() != NULL); assertion(_cyclicCommRight->isConnected());

		_cyclicCommLeft = leftComm;
		_cyclicCommRight = rightComm;
	}
}

void ParallelMatrixOperations::multiply
(
	TarchMatrix& leftMatrix,
	TarchMatrix& rightMatrix,
	TarchMatrix& result,
	std::vector<int>& offsets,
	int p, int q, int r)
{
	assertion2(result.rows() == leftMatrix.rows(), result.rows(), leftMatrix.rows());
	assertion2(result.cols() == rightMatrix.cols(), result.cols(), rightMatrix.cols());
	assertion2(leftMatrix.cols() == rightMatrix.rows(), leftMatrix.cols(), rightMatrix.rows());

	// if serial computation on single processor, i.e, no master-slave mode
	if( not utils::MasterSlave::_masterMode && not utils::MasterSlave::_slaveMode){
		multiply(leftMatrix, rightMatrix, result);

	// if parallel computation on p processors, i.e., master-slave mode
	}else{
		assertion(utils::MasterSlave::_communication.get() != NULL);
		assertion(utils::MasterSlave::_communication->isConnected());

		// The result matrix is of size (p x r)
		// if p equals r (and p = global_n), we have to perform the
		// cyclic communication with block-wise matrix-matrix multiplication
		if(p == r){
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
	EigenMatrix& leftMatrix,
	EigenMatrix& rightMatrix,
	EigenMatrix& result,
	std::vector<int>& offsets,
	int p, int q, int r)
{
	assertion2(result.rows() == leftMatrix.rows(), result.rows(), leftMatrix.rows());
	assertion2(result.cols() == rightMatrix.cols(), result.cols(), rightMatrix.cols());
	assertion2(leftMatrix.cols() == rightMatrix.rows(), leftMatrix.cols(), rightMatrix.rows());

	// if serial computation on single processor, i.e, no master-slave mode
	if( not utils::MasterSlave::_masterMode && not utils::MasterSlave::_slaveMode){
		result = leftMatrix * rightMatrix;

	// if parallel computation on p processors, i.e., master-slave mode
	}else{
		assertion(utils::MasterSlave::_communication.get() != NULL);
		assertion(utils::MasterSlave::_communication->isConnected());

		// The result matrix is of size (p x r)
		// if p equals r (and p = global_n), we have to perform the
		// cyclic communication with block-wise matrix-matrix multiplication
		if(p == r){
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

// @brief multiplies matrices based on a dot-product computation with a rectangular result matrix
void ParallelMatrixOperations::_multiplyNM(
		TarchMatrix& leftMatrix, TarchMatrix& rightMatrix, TarchMatrix& result, std::vector<int>& offsets, int p, int q, int r)
{

  for(int i = 0; i < leftMatrix.rows(); i++){
	  int rank = 0;
	  // find rank of processor that stores the result
	  // the second while is necessary if processors with no vertices are present
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

		  // TODO: better: implement a reduce-operation (no loop over all slaves)
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

// @brief multiplies matrices based on a dot-product computation with a rectangular result matrix
void ParallelMatrixOperations::_multiplyNM(
		EigenMatrix& leftMatrix, EigenMatrix& rightMatrix, EigenMatrix& result, std::vector<int>& offsets, int p, int q, int r)
{

}

// @brief multiplies matrices based on a cyclic communication and block-wise matrix multiplication with a quadratic result matrix
void ParallelMatrixOperations::_multiplyNN(
		TarchMatrix& leftMatrix, TarchMatrix& rightMatrix, TarchMatrix& result, std::vector<int>& offsets, int p, int q, int r)
{

}

// @brief multiplies matrices based on a cyclic communication and block-wise matrix multiplication with a quadratic result matrix
void ParallelMatrixOperations::_multiplyNN(
		EigenMatrix& leftMatrix, EigenMatrix& rightMatrix, EigenMatrix& result, std::vector<int>& offsets, int p, int q, int r)
{

}


void ParallelMatrixOperations::multiply
(
	TarchMatrix& A,
	TarchVector& v,
	TarchVector& result,
	std::vector<int>& offsets,
	int p, int q)
{

}

void ParallelMatrixOperations::multiply
(
	EigenMatrix& A,
	EigenVector& v,
	EigenVector& result,
	std::vector<int>& offsets,
	int p, int q)
{

}


}}} // namespace precice, cplscheme, impl

#endif
