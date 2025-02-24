#pragma once

#ifndef PRECICE_NO_MPI

#include <Eigen/Core>
#include <memory>
#include <stddef.h>
#include <string>
#include <vector>

#include "com/Communication.hpp"
#include "com/MPIPortsCommunication.hpp"
#include "com/Request.hpp"
#include "com/SharedPointer.hpp"
#include "logging/LogMacros.hpp"
#include "logging/Logger.hpp"
#include "precice/impl/Types.hpp"
#include "utils/IntraComm.hpp"
#include "utils/assertion.hpp"

namespace precice::acceleration::impl {

class ParallelMatrixOperations {
public:
  ~ParallelMatrixOperations();

  /// Initializes the acceleration.
  void initialize(const bool needCyclicComm);

  template <typename Derived1, typename Derived2>
  void multiply(
      Eigen::PlainObjectBase<Derived1> &leftMatrix,
      Eigen::PlainObjectBase<Derived2> &rightMatrix,
      Eigen::PlainObjectBase<Derived2> &result,
      const std::vector<int> &          offsets,
      int p, int q, int r,
      bool cyclicComm            = true,
      bool dotProductComputation = true)
  {
    PRECICE_TRACE();
    PRECICE_ASSERT(result.cols() == rightMatrix.cols(), result.cols(), rightMatrix.cols());
    PRECICE_ASSERT(leftMatrix.cols() == rightMatrix.rows(), leftMatrix.cols(), rightMatrix.rows());

    // if serial computation on single processor
    if (!utils::IntraComm::isParallel()) {
      result.noalias() = leftMatrix * rightMatrix;

      // if parallel computation on p processors
    } else {
      PRECICE_ASSERT(utils::IntraComm::getCommunication() != nullptr);
      PRECICE_ASSERT(utils::IntraComm::getCommunication()->isConnected());

      // The result matrix is of size (p x r)
      // if the cyclic communication is needed, we use block-wise matrix-matrix multiplication
      if (cyclicComm) {
        PRECICE_ASSERT(_needCyclicComm);
        PRECICE_ASSERT(_cyclicCommLeft.get() != nullptr);
        PRECICE_ASSERT(_cyclicCommLeft->isConnected());
        PRECICE_ASSERT(_cyclicCommRight.get() != nullptr);
        PRECICE_ASSERT(_cyclicCommRight->isConnected());

        _multiply_cyclic(leftMatrix, rightMatrix, result, offsets, p, q, r);

        // case the cyclic communication is not needed, i.e., usually p = number of columns of the least squares system
        // perform parallel multiplication based on dot-product
      } else {
        if (dotProductComputation)
          _multiplyNM_dotProduct(leftMatrix, rightMatrix, result, offsets, p, q, r);
        else
          _multiplyNM_block(leftMatrix, rightMatrix, result, offsets, p, q, r);
      }
    }
  }

  /** @brief: Method computes the matrix-matrix/matrix-vector product of a (p x q)
    * matrix that is distributed column-wise (e.g. pseudoInverse Z), with a matrix/vector
    * of size (q x r) with r=1 or r=cols, that is distributed row-wise (e.g. _matrixW, _matrixV, residual).
    *
    * In each case mat-mat or mat-vec product, the result is of size (m x m) or (m x 1), where
    * m is the number of cols, i.e., small such that the result is stored on each proc.
    *
    * @param[in] p - first dimension, i.e., overall (global) number cols of result matrix
    * @param[in] q - inner dimension
    * @param[in] r - second dimension, i.e., 1 or overall (global) number cols of result matrix
    *
    */
  template <typename Derived1, typename Derived2, typename Derived3>
  void multiply(
      const Eigen::MatrixBase<Derived1> &leftMatrix,
      const Eigen::MatrixBase<Derived2> &rightMatrix,
      Eigen::PlainObjectBase<Derived3> & result,
      int p, int q, int r)
  {
    PRECICE_TRACE();
    PRECICE_ASSERT(leftMatrix.rows() == p, leftMatrix.rows(), p);
    PRECICE_ASSERT(leftMatrix.cols() == rightMatrix.rows(), leftMatrix.cols(), rightMatrix.rows());
    PRECICE_ASSERT(result.rows() == p, result.rows(), p);
    PRECICE_ASSERT(result.cols() == r, result.cols(), r);

    Eigen::MatrixXd localResult(result.rows(), result.cols());
    localResult.noalias() = leftMatrix * rightMatrix;

    // if serial computation on single processor
    if (!utils::IntraComm::isParallel()) {
      result = localResult;
    } else {
      utils::IntraComm::allreduceSum(localResult, result);
    }
  }

private:
  logging::Logger _log{"acceleration::ParallelMatrixOperations"};

  // @brief multiplies matrices based on a cyclic communication and block-wise matrix multiplication
  template <typename Derived1, typename Derived2>
  void _multiply_cyclic(
      Eigen::PlainObjectBase<Derived1> &leftMatrix,
      Eigen::PlainObjectBase<Derived2> &rightMatrix,
      Eigen::PlainObjectBase<Derived2> &result,
      const std::vector<int> &          offsets,
      int p, int q, int r)
  {
    PRECICE_TRACE();
    /*
     * For multiplication W_til * Z = J
     * -----------------------------------------------------------------------
     * p = n_global, q = m, r = n_global_primary
     *
     * leftMatrix:  local: (n_local x m) 		global: (n_global x m)
     * rightMatrix: local: (m x n_local_primary) 		global: (m x n_global_primary)
     * result: 		local: (n_global x n_local_primary) global: (n_global x n_global_primary)
     * -----------------------------------------------------------------------
     */

    PRECICE_ASSERT(_needCyclicComm);
    PRECICE_ASSERT(leftMatrix.cols() == q, leftMatrix.cols(), q);
    PRECICE_ASSERT(result.rows() == p, result.rows(), p);

    //int nextProc = (utils::IntraComm::getRank() + 1) % utils::IntraComm::getSize();
    int prevProc = (utils::IntraComm::getRank() - 1 < 0) ? utils::IntraComm::getSize() - 1 : utils::IntraComm::getRank() - 1;
    int rows_rcv = (prevProc > 0) ? offsets[prevProc + 1] - offsets[prevProc] : offsets[1];
    //Eigen::MatrixXd leftMatrix_rcv = Eigen::MatrixXd::Zero(rows_rcv, q);
    Eigen::MatrixXd leftMatrix_rcv(rows_rcv, q);

    com::PtrRequest requestSend;
    com::PtrRequest requestRcv;

    // initiate asynchronous send operation of leftMatrix (W_til) --> nextProc (this data is needed in cycle 1)    dim: n_local x cols
    if (leftMatrix.size() > 0)
      requestSend = _cyclicCommRight->aSend(leftMatrix, 0);

    // initiate asynchronous receive operation for leftMatrix (W_til) from previous processor --> W_til      dim: rows_rcv x cols
    if (leftMatrix_rcv.size() > 0)
      requestRcv = _cyclicCommLeft->aReceive(leftMatrix_rcv, 0);

    // compute diagonal blocks where all data is local and no communication is needed
    // compute block matrices of J_inv of size (n_til x n_til), n_til = local n
    Eigen::MatrixXd diagBlock(leftMatrix.rows(), leftMatrix.rows());
    diagBlock.noalias() = leftMatrix * rightMatrix;

    // set block at corresponding row-index on proc
    int off = offsets[utils::IntraComm::getRank()];
    PRECICE_ASSERT(result.cols() == diagBlock.cols(), result.cols(), diagBlock.cols());
    result.block(off, 0, diagBlock.rows(), diagBlock.cols()) = diagBlock;

    /**
		 * cyclic send-receive operation
		 */
    for (int cycle = 1; cycle < utils::IntraComm::getSize(); cycle++) {

      // wait until W_til from previous processor is fully received
      if (requestSend != nullptr)
        requestSend->wait();
      if (requestRcv != nullptr)
        requestRcv->wait();

      // leftMatrix (leftMatrix_rcv) is available - needed for local multiplication and hand over to next proc
      Eigen::MatrixXd leftMatrix_copy(leftMatrix_rcv);

      // initiate async send to hand over leftMatrix (W_til) to the next proc (this data will be needed in the next cycle)    dim: n_local x cols
      if (cycle < utils::IntraComm::getSize() - 1) {
        if (leftMatrix_copy.size() > 0)
          requestSend = _cyclicCommRight->aSend(leftMatrix_copy, 0);
      }

      // compute proc that owned leftMatrix_rcv (Wtil_rcv) at the very beginning for each cycle
      int sourceProc_nextCycle = (utils::IntraComm::getRank() - (cycle + 1) < 0) ? utils::IntraComm::getSize() + (utils::IntraComm::getRank() - (cycle + 1)) : utils::IntraComm::getRank() - (cycle + 1);

      int sourceProc = (utils::IntraComm::getRank() - cycle < 0) ? utils::IntraComm::getSize() + (utils::IntraComm::getRank() - cycle) : utils::IntraComm::getRank() - cycle;

      int rows_rcv_nextCycle = (sourceProc_nextCycle > 0) ? offsets[sourceProc_nextCycle + 1] - offsets[sourceProc_nextCycle] : offsets[1];
      rows_rcv               = (sourceProc > 0) ? offsets[sourceProc + 1] - offsets[sourceProc] : offsets[1];
      leftMatrix_rcv         = Eigen::MatrixXd::Zero(rows_rcv_nextCycle, q);

      // initiate asynchronous receive operation for leftMatrix (W_til) from previous processor --> W_til (this data is needed in the next cycle)
      if (cycle < utils::IntraComm::getSize() - 1) {
        if (leftMatrix_rcv.size() > 0) // only receive data, if data has been sent
          requestRcv = _cyclicCommLeft->aReceive(leftMatrix_rcv, 0);
      }

      if (requestSend != nullptr)
        requestSend->wait();
      // compute block with new local data
      Eigen::MatrixXd block(rows_rcv, rightMatrix.cols());
      block.noalias() = leftMatrix_copy * rightMatrix;

      // set block at corresponding index in J_inv
      // the row-offset of the current block is determined by the proc that sends the part of the W_til matrix
      // note: the direction and ordering of the cyclic sending operation is chosen s.t. the computed block is
      //       local on the current processor (in J_inv).
      off = offsets[sourceProc];
      PRECICE_ASSERT(result.cols() == block.cols(), result.cols(), block.cols());
      result.block(off, 0, block.rows(), block.cols()) = block;
    }
  }

  // @brief multiplies matrices based on a dot-product computation with a rectangular result matrix
  template <typename Derived1, typename Derived2>
  void _multiplyNM_dotProduct(
      Eigen::PlainObjectBase<Derived1> &leftMatrix,
      Eigen::PlainObjectBase<Derived2> &rightMatrix,
      Eigen::PlainObjectBase<Derived2> &result,
      const std::vector<int> &          offsets,
      int p, int q, int r)
  {
    PRECICE_TRACE();
    for (int i = 0; i < leftMatrix.rows(); i++) {
      Rank rank = 0;
      // find rank of processor that stores the result
      // the second while is necessary if processors with no vertices are present
      // Note: the >'=' here is crucial: In case some procs do not have any vertices,
      // this while loop continues incrementing rank if entries in offsets are equal, i.e.,
      // it runs to the next non-empty proc.
      while (i >= offsets[rank + 1])
        rank++;

      Eigen::VectorXd lMRow = leftMatrix.row(i);

      for (int j = 0; j < r; j++) {

        Eigen::VectorXd rMCol  = rightMatrix.col(j);
        double          res_ij = utils::IntraComm::dot(lMRow, rMCol);

        // find proc that needs to store the result.
        int local_row;
        if (utils::IntraComm::getRank() == rank) {
          local_row            = i - offsets[rank];
          result(local_row, j) = res_ij;
        }
      }
    }
  }

  /// Multiplies matrices based on a SAXPY-like block-wise computation with a rectangular result matrix of dimension n x m
  template <typename Derived1, typename Derived2>
  void _multiplyNM_block(
      Eigen::PlainObjectBase<Derived1> &leftMatrix,
      Eigen::PlainObjectBase<Derived2> &rightMatrix,
      Eigen::PlainObjectBase<Derived2> &result,
      const std::vector<int> &          offsets,
      int p, int q, int r)
  {
    PRECICE_TRACE();

    // ensure that both matrices are stored in the same order. Important for reduce function, that adds serialized data.
    PRECICE_ASSERT(static_cast<int>(leftMatrix.IsRowMajor) == static_cast<int>(rightMatrix.IsRowMajor),
                   static_cast<int>(leftMatrix.IsRowMajor), static_cast<int>(rightMatrix.IsRowMajor));

    // multiply local block (saxpy-based approach)
    // dimension: (n_global x n_local) * (n_local x m) = (n_global x m)
    Eigen::MatrixXd block = Eigen::MatrixXd::Zero(p, r);
    block.noalias()       = leftMatrix * rightMatrix;

    // all blocks have size (n_global x m)
    // Note: if procs have no vertices, the block size remains (n_global x m), however,
    // 	     it must be initialized with zeros, so zeros are added for those procs)

    // sum up blocks on the primary rank, reduce
    Eigen::MatrixXd summarizedBlocks = Eigen::MatrixXd::Zero(p, r); /// @todo: only primary rank should allocate memory.
    utils::IntraComm::reduceSum(block, summarizedBlocks);

    // secondary ranks wait to receive their local result
    if (utils::IntraComm::isSecondary()) {
      if (result.size() > 0)
        utils::IntraComm::getCommunication()->receive(result, 0);
    }

    // primary rank distributes the sub blocks of the results
    if (utils::IntraComm::isPrimary()) {
      // distribute blocks of summarizedBlocks (result of multiplication) to corresponding secondary ranks
      result = summarizedBlocks.block(0, 0, offsets[1], r);

      for (Rank secondaryRank : utils::IntraComm::allSecondaryRanks()) {
        int off       = offsets[secondaryRank];
        int send_rows = offsets[secondaryRank + 1] - offsets[secondaryRank];

        if (summarizedBlocks.block(off, 0, send_rows, r).size() > 0) {
          // necessary to save the matrix-block that is to be sent in a temporary matrix-object
          // otherwise, the send routine walks over the bounds of the block (matrix structure is still from the entire matrix)
          Eigen::MatrixXd sendBlock = summarizedBlocks.block(off, 0, send_rows, r);
          utils::IntraComm::getCommunication()->send(sendBlock, secondaryRank);
        }
      }
    }
  }

  /// Communication between neighboring ranks, backwards
  com::PtrCommunication _cyclicCommLeft = nullptr;

  /// Communication between neighboring ranks, forward
  com::PtrCommunication _cyclicCommRight = nullptr;

  bool _needCyclicComm = true;

  /** Establishes the circular connection between ranks
   *
   * This creates and connects the ranks.
   *
   * @precondition _cyclicCommLeft and _cyclicCommRight must be nullptr
   * @postcondition _cyclicCommLeft, _cyclicCommRight are connected
   */
  void establishCircularCommunication();

  /** Closes the circular connection between ranks
   *
   * @precondition establishCircularCommunication() was called
   * @postcondition _cyclicCommLeft, _cyclicCommRight are disconnected and set to nullptr
   */
  void closeCircularCommunication();
};

} // namespace precice::acceleration::impl

#endif
