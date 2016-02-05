/*
 * SVDFactorization.hpp
 *
 *  Created on: Feb 5, 2016
 *      Author: scheufks
 */

#ifndef SVDFACTORIZATION_HPP_
#define SVDFACTORIZATION_HPP_

#include "mesh/SharedPointer.hpp"
#include "tarch/logging/Log.h"
#include "utils/MasterSlave.hpp"
#include <Eigen/Dense>
#include <limits>
#include <deque>
#include <fstream>



// ------- CLASS DEFINITION

namespace precice {
namespace cplscheme {
namespace impl {

/**
 * @brief Class that provides functionality to maintain a SVD decomposition of a matrix
 * via succesive rank-1 updates and truncation with respect to the truncation threshold eps.
 */
class SVDFactorization
{
public:
  // Eigen
  typedef Eigen::MatrixXd Matrix;
  typedef Eigen::VectorXd Vector;



  /**
   * @brief Constructor.
   */
  SVDFactorization(
      double eps,
      PtrPreconditioner preconditioner);

   /**
    * @brief Destructor, empty.
    */
   virtual ~SVDFactorization() {}

   /** @brief: updates the SVD decomposition with the rank-1 update A*B^T, i.e.,
    *               _psi * _sigma * _phi^T + A*B^T
    *  and overrides the internal SVD representation. After the update, the SVD is
    *  truncated according to the threshold _truncationEps
    *
    *  THIS METHOD RUNS IN SERIAL ON MASTER PROC!
    */
   void update(Matrix& A, Matrix& B);


   /**
    * @brief resets the SVD factorization
    */
   void reset();

   /**
    * @brief returns a matrix representation of the orthogonal matrix Psi, A = Psi * Sigma * Phi^T
    */
   Matrix& matrixPsi();

   /**
    * @brief returns a matrix representation of the orthogonal matrix Sigma, A = Psi * Sigma * Phi^T
    */
   Matrix& matrixSigma();

   /**
    * @brief returns a matrix representation of the orthogonal matrix Phi, A = Psi * Sigma * Phi^T
    */
   Matrix& matrixPhi();

   // @brief returns the number of columns in the QR-decomposition
   int cols();
   // @brief returns the number of rows in the QR-decomposition
   int rows();

   void setPrecondApplied(bool b);

   bool isPrecondApplied();

   // @brief optional file-stream for logging output
   void setfstream(std::fstream* stream);

private:


  // @brief Logging device.
  static tarch::logging::Log _log;

  /// @brief preconditioner for least-squares system if vectorial system is used.
  PtrPreconditioner _preconditioner;

  /// @brief: SVD factorization of the matrix J = _psi * _sigma * _phi^T
  Matrix _psi;
  Matrix _phi;
  Matrix _sigma;

  int _rows;
  int _cols;

  ///@brief: Truncation parameter for the updated SVD decomposition
  double _truncationEps;

  /// @brief: true if the preconditioner has been applied appropriate to the updated SVD decomposition
  bool   _preconditionerApplied;

  // @brief optional infostream that writes information to file
  std::fstream* _infostream;
  bool _fstream_set;

};

}}} // namespace precice, cplscheme, impl



#endif /* SVDFACTORIZATION_HPP_ */
