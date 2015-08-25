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
   void multiply(EigenMatrix& leftMatrix,
		   	   	 EigenMatrix& rightMatrix,
		   	   	 EigenMatrix& result,
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
	* @brief multiplies Eigen matrix with Eigen vector in parallel or serial execution.
	* 		 The multiplication is based on a dot-product computation in parallel execution.
	* @param [IN] p - first dimension, i.e., overall (global) number of rows
	* @param [IN] q - second dimension, i.e., overall (global) number cols
	*/
/*   void multiply(EigenMatrix& A,
		   	     EigenVector& v,
		   	     EigenVector& result,
		   	     std::vector<int>& offsets,
		   	     int p, int q);
*/

private:

   // @brief Logging device.
   static tarch::logging::Log _log;

   // @brief multiplies matrices based on a dot-product computation with a rectangular result matrix
   void _multiplyNM(TarchMatrix& leftMatrix, TarchMatrix& rightMatrix, TarchMatrix& result, std::vector<int>& offsets, int p, int q, int r);
   // @brief multiplies matrices based on a dot-product computation with a rectangular result matrix
   void _multiplyNM(EigenMatrix& leftMatrix, EigenMatrix& rightMatrix, EigenMatrix& result, std::vector<int>& offsets, int p, int q, int r);

   // @brief multiplies matrices based on a cyclic communication and block-wise matrix multiplication with a quadratic result matrix
   void _multiplyNN(TarchMatrix& leftMatrix, TarchMatrix& rightMatrix, TarchMatrix& result, std::vector<int>& offsets, int p, int q, int r);
   // @brief multiplies matrices based on a cyclic communication and block-wise matrix multiplication with a quadratic result matrix
   void _multiplyNN(EigenMatrix& leftMatrix, EigenMatrix& rightMatrix, EigenMatrix& result, std::vector<int>& offsets, int p, int q, int r);

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
