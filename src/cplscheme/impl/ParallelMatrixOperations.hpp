/*
 * ParallelMatrixOperations.hpp
 *
 *  Created on: Aug 21, 2015
 *      Author: scheufks
 */
#ifndef PRECICE_NO_MPI
#ifndef PARALLELMATRIXOPERATIONS_HPP_
#define PARALLELMATRIXOPERATIONS_HPP_


#include "tarch/logging/Log.h"
#include "tarch/la/DynamicColumnMatrix.h"
#include "tarch/la/DynamicMatrix.h"
#include "tarch/la/DynamicVector.h"
#include "com/MPIPortsCommunication.hpp"
#include "com/Communication.hpp"
#include "utils/MasterSlave.hpp"
#include "utils/Globals.hpp"
#include "Eigen/Dense"

namespace precice {
namespace cplscheme {
namespace impl {

class ParallelMatrixOperations
{
public:

	// tarch
	typedef tarch::la::DynamicVector<double> TarchVector;
	typedef tarch::la::DynamicMatrix<double> TarchMatrix;
	typedef tarch::la::DynamicColumnMatrix<double> TarchColumnMatrix;

	// Eigen
	typedef Eigen::MatrixXd EigenMatrix;
	typedef Eigen::VectorXd EigenVector;

  /**
   * @brief Constructor.
   */
	ParallelMatrixOperations ();

   /**
    * @brief Destructor, empty.
    */
   virtual ~ParallelMatrixOperations(){};


   /**
    * @brief Initializes the post-processing.
    */
   void initialize(com::Communication::SharedPointer leftComm,
		   	   	   com::Communication::SharedPointer rightComm);

   /**
    * @brief multiplies tarch matrices in parallel or serial execution.
    * 		 This class is specialized for multiplication of matrices in the
    * 		 coupling scheme.
    * 		 If leftMatrix.rows() = n_global and the result matrix is of size
    * 		 (n_global x n_global) a cyclic communication is used and the overall
    * 		 matrix is computed block-wise and distributed.
    * 		 If the result matrix is of size (n_global x LS_cols) the multiplication
    * 		 is based on a dot-product computation.
    * @param [IN] p - first dimension, i.e., overall (global) number of rows
    * @param [IN] q - inner dimension
    * @param [IN] r - second dimension, i.e., overall (global) number cols of result matrix
    */
   void multiply(TarchMatrix& leftMatrix,
		   	   	 TarchMatrix& rightMatrix,
		   	   	 TarchMatrix& result,
		   	   	 std::vector<int>& offsets,
		   	   	 int p, int q, int r);
   void multiply(TarchMatrix& leftMatrix,
   		   	   	 TarchColumnMatrix& rightMatrix,
   		   	   	 TarchMatrix& result,
   		   	   	 std::vector<int>& offsets,
   		   	   	 int p, int q, int r);
   /**
    * @brief multiplies tarch matrix with tarch vector in parallel or serial execution.
    * 		 The multiplication is based on a dot-product computation in parallel execution.
    * @param [IN] p - first dimension, i.e., overall (global) number of rows
    * @param [IN] q - second dimension, i.e., overall (global) number cols
    */
   void multiply(TarchMatrix& A,
		   	   	 TarchVector& v,
		   	   	 TarchVector& result,
		   	   	 std::vector<int>& offsets,
		   	   	 int p, int q);

   /**
	* @brief multiplies Eigen matrices in parallel or serial execution.
	* 		 This class is specialized for multiplication of matrices in the
	* 		 coupling scheme.
	* 		 If leftMatrix.rows() = n_global and the result matrix is of size
	* 		 (n_global x n_global) a cyclic communication is used and the overall
	* 		 matrix is computed block-wise and distributed.
	* 		 If the result matrix is of size (n_global x LS_cols) the multiplication
	* 		 is based on a dot-product computation.
	* @param [IN] p - first dimension, i.e., overall (global) number of rows
    * @param [IN] q - inner dimension
    * @param [IN] r - second dimension, i.e., overall (global) number cols of result matrix
	*/
   template<typename Derived1, typename Derived2>
   void multiply(Eigen::PlainObjectBase<Derived1>& leftMatrix,
		   	   	 Eigen::PlainObjectBase<Derived2>& rightMatrix,
		   	   	 Eigen::PlainObjectBase<Derived2>& result,
   		   	   	 std::vector<int>& offsets,
   		   	   	 int p, int q, int r)
   {
		preciceTrace("multiply()");
		assertion2(result.cols() == rightMatrix.cols(), result.cols(), rightMatrix.cols());
		assertion2(leftMatrix.cols() == rightMatrix.rows(), leftMatrix.cols(), rightMatrix.rows());

		// if serial computation on single processor, i.e, no master-slave mode
		if( not utils::MasterSlave::_masterMode && not utils::MasterSlave::_slaveMode){
			result.noalias() = leftMatrix * rightMatrix;

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


private:

   // @brief Logging device.
   static tarch::logging::Log _log;

   // @brief multiplies matrices based on a dot-product computation with a rectangular result matrix
   void _multiplyNM(TarchMatrix& leftMatrix, TarchMatrix& rightMatrix, TarchMatrix& result, std::vector<int>& offsets, int p, int q, int r);

   // @brief multiplies matrices based on a cyclic communication and block-wise matrix multiplication with a quadratic result matrix
   void _multiplyNN(TarchMatrix& leftMatrix, TarchMatrix& rightMatrix, TarchMatrix& result, std::vector<int>& offsets, int p, int q, int r);



   // @brief multiplies matrices based on a cyclic communication and block-wise matrix multiplication with a quadratic result matrix
   template<typename Derived1, typename Derived2>
   void _multiplyNN(Eigen::PlainObjectBase<Derived1>& leftMatrix,
		   	   	    Eigen::PlainObjectBase<Derived2>& rightMatrix,
		   	   	    Eigen::PlainObjectBase<Derived2>& result,
		   	   	    std::vector<int>& offsets,
		   	   	    int p, int q, int r)
   {
		preciceTrace("multiplyNN()");
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

		assertion2(leftMatrix.cols() == q, leftMatrix.cols(), q);
		assertion2(leftMatrix.rows() == rightMatrix.cols(), leftMatrix.rows(), rightMatrix.cols());
		assertion2(result.rows() == p, result.rows(), p);

		//int nextProc = (utils::MasterSlave::_rank + 1) % utils::MasterSlave::_size;
		int prevProc = (utils::MasterSlave::_rank -1 < 0) ? utils::MasterSlave::_size-1 : utils::MasterSlave::_rank -1;
		int rows_rcv = (prevProc > 0) ? offsets[prevProc+1] - offsets[prevProc] : offsets[1];
		//EigenMatrix leftMatrix_rcv = EigenMatrix::Zero(rows_rcv, q);
		EigenMatrix leftMatrix_rcv(rows_rcv, q);

		com::Request::SharedPointer requestSend;
		com::Request::SharedPointer requestRcv;

		// initiate asynchronous send operation of leftMatrix (W_til) --> nextProc (this data is needed in cycle 1)    dim: n_local x cols
		if(leftMatrix.size() > 0)
			requestSend = _cyclicCommRight->aSend(leftMatrix.data(), leftMatrix.size(), 0);

		// initiate asynchronous receive operation for leftMatrix (W_til) from previous processor --> W_til      dim: rows_rcv x cols
		if(leftMatrix_rcv.size() > 0)
			requestRcv = _cyclicCommLeft->aReceive(leftMatrix_rcv.data(), leftMatrix_rcv.size(), 0);

		// compute diagonal blocks where all data is local and no communication is needed
		// compute block matrices of J_inv of size (n_til x n_til), n_til = local n
		EigenMatrix diagBlock(leftMatrix.rows(), leftMatrix.rows());
		diagBlock.noalias() = leftMatrix * rightMatrix;

		// set block at corresponding row-index on proc
		int off = offsets[utils::MasterSlave::_rank];
		assertion2(result.cols() == diagBlock.cols(), result.cols(), diagBlock.cols());
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
			EigenMatrix leftMatrix_copy(leftMatrix_rcv);

			// initiate async send to hand over leftMatrix (W_til) to the next proc (this data will be needed in the next cycle)    dim: n_local x cols
			if(cycle < utils::MasterSlave::_size-1){
			  if(leftMatrix_copy.size() > 0)
				  requestSend = _cyclicCommRight->aSend(leftMatrix_copy.data(), leftMatrix_copy.size(), 0);
			}

			// compute proc that owned leftMatrix_rcv (Wtil_rcv) at the very beginning for each cylce
			int sourceProc_nextCycle = (utils::MasterSlave::_rank - (cycle+1) < 0) ?
				  utils::MasterSlave::_size + (utils::MasterSlave::_rank - (cycle+1)) : utils::MasterSlave::_rank - (cycle+1);

			int sourceProc = (utils::MasterSlave::_rank - cycle < 0) ?
				  utils::MasterSlave::_size + (utils::MasterSlave::_rank - cycle) : utils::MasterSlave::_rank - cycle;

			int rows_rcv_nextCycle = (sourceProc_nextCycle > 0) ? offsets[sourceProc_nextCycle+1] - offsets[sourceProc_nextCycle] : offsets[1];
			rows_rcv = (sourceProc > 0) ? offsets[sourceProc+1] - offsets[sourceProc] : offsets[1];
			leftMatrix_rcv = EigenMatrix::Zero(rows_rcv_nextCycle, q);


			// initiate asynchronous receive operation for leftMatrix (W_til) from previous processor --> W_til (this data is needed in the next cycle)
			if(cycle < utils::MasterSlave::_size-1){
			  if(leftMatrix_rcv.size() > 0) // only receive data, if data has been sent
				  requestRcv = _cyclicCommLeft->aReceive(leftMatrix_rcv.data(), leftMatrix_rcv.size(), 0);
			}

			// compute block with new local data
			EigenMatrix block(rows_rcv, rightMatrix.cols());
			block.noalias() = leftMatrix_copy * rightMatrix;

			// set block at corresponding index in J_inv
			// the row-offset of the current block is determined by the proc that sends the part of the W_til matrix
			// note: the direction and ordering of the cyclic sending operation is chosen s.t. the computed block is
			//       local on the current processor (in J_inv).
			off = offsets[sourceProc];
			assertion2(result.cols() == block.cols(), result.cols(), block.cols());
			for(int ii = 0; ii < block.rows(); ii++)
			  for(int jj = 0; jj < result.cols(); jj++)
			  {
				  result(ii+off, jj) = block(ii, jj);
			  }
		}
   }

   // @brief multiplies matrices based on a dot-product computation with a rectangular result matrix
   template<typename Derived1, typename Derived2>
   void _multiplyNM(Eigen::PlainObjectBase<Derived1>& leftMatrix,
		   	   	    Eigen::PlainObjectBase<Derived2>& rightMatrix,
		   	   	    Eigen::PlainObjectBase<Derived2>& result,
		   	   	    std::vector<int>& offsets,
		   	   	    int p, int q, int r)
   {
		preciceTrace("multiplyNM()");
		for(int i = 0; i < leftMatrix.rows(); i++){
		  int rank = 0;
		  // find rank of processor that stores the result
		  // the second while is necessary if processors with no vertices are present
		  // Note: the >'=' here is crucial: In case some procs do not have any vertices,
		  // this while loop continues incrementing rank if entries in offsets are equal, i.e.,
		  // it runs to the next non-empty proc.
		  while(i >= offsets[rank+1]) rank++;

		  EigenVector lMRow = leftMatrix.row(i);

		  for(int j = 0; j < r; j++){

			  EigenVector rMCol = rightMatrix.col(j);
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


	/**
	* @brief Communication between neighboring slaves, backwards
	*/
	com::Communication::SharedPointer _cyclicCommLeft;

	/**
	* @brief Communication between neighboring slaves, forward
	*/
	com::Communication::SharedPointer _cyclicCommRight;

};

}}} // namespace precice, cplscheme, impl


#endif /* PARALLELMATRIXOPERATIONS_HPP_ */
#endif
