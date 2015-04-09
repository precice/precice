// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_QRFACTORIZATION_HPP_
#define PRECICE_QRFACTORIZATION_HPP_

#include "mesh/SharedPointer.hpp"
#include "tarch/logging/Log.h"
#include "tarch/la/DynamicColumnMatrix.h"
#include "tarch/la/DynamicMatrix.h"
#include "tarch/la/DynamicVector.h"
#include <Eigen/Dense>
#include <limits>
#include <deque>
#include <fstream>



// ------- CLASS DEFINITION

namespace precice {
namespace cplscheme {
namespace impl {

/**
 * @brief Base Class for quasi-Newton post processing schemes
 * 
 */
class QRFactorization
{
public:
  
  // tarch
  typedef tarch::la::DynamicVector<double> DataValues;
  typedef tarch::la::DynamicColumnMatrix<double> DataMatrix;
  typedef tarch::la::DynamicMatrix<double> Matrix;
   
  // Eigen
  typedef Eigen::MatrixXd EigenMatrix;
  typedef Eigen::VectorXd EigenVector;

  /**
   * @brief Constructor.
   */
   QRFactorization (
      DataMatrix A,
      double omega=0,
      double theta=1./0.7,
      double sigma=std::numeric_limits<double>::min()
 		    );
   
   /**
   * @brief Constructor.
   */
   QRFactorization (
      EigenMatrix A,
      double omega=0,
      double theta=1./0.7,
      double sigma=std::numeric_limits<double>::min()
 		    );
   
    /**
   * @brief Constructor.
   */
   QRFactorization (
      EigenMatrix Q,
      EigenMatrix R,
      size_t rows,
      size_t cols,
      double omega=0,
      double theta=1./0.7,
      double sigma=std::numeric_limits<double>::min()
 		    );
   
   /**
    * @brief Destructor, empty.
    */
   virtual ~QRFactorization() {}


   /**
    * @brief resets the QR factorization zo zero Q(0:0, 0:0)R(0:0, 0:0)
    */
   void reset();
   
   /**
    * @brief resets the QR factorization to the given factorization Q, R
    */
   void reset(
	EigenMatrix Q, 
	EigenMatrix R, 
	size_t rows, 
	size_t cols, 
	double omega=0,
        double theta=1./0.7,
        double sigma=std::numeric_limits<double>::min());
   
   /**
    * @brief resets the QR factorization to be the factorization of A = QR
    */
   void reset(
	EigenMatrix A,
	double omega=0,
        double theta=1./0.7,
        double sigma=std::numeric_limits<double>::min());
   
   void insertColumn(size_t k, EigenVector& v);
   void insertColumn(size_t k, DataValues& v);
   
   /**
   * updates the factorization A=Q[1:n,1:m]R[1:m,1:n] when the kth column of A is deleted. 
   * Returns the deleted column v(1:n)
   */
   void deleteColumn(size_t k);
   
   EigenMatrix& matrixQ();
   EigenMatrix& matrixR();
   
   size_t cols();
   size_t rows();

private:
  
  struct givensRot{
    size_t i, j;
    double sigma, gamma;
  };
   
  /**
  * @short assuming Q(1:n,1:m) has nearly orthonormal columns, this procedure
  *   orthogonlizes v(1:n) to the columns of Q, and normalizes the result.
  *   r(1:n) is the array of Fourier coefficients, and rho is the distance
  *   from v to range of Q, r and its corrections are computed in double
  *   precision.
  */
  int orthogonalize(EigenVector& v, EigenVector& r, double &rho, size_t colNum);
  
  /**
  * @short computes parameters for givens matrix G for which  (x,y)G = (z,0). replaces (x,y) by (z,0)
  */
  void computeReflector(givensRot &grot, double &x, double &y);
  
  /**
  *  @short this procedure replaces the two column matrix [p(k:l-1), q(k:l-1)] by [p(k:l), q(k:l)]*G, 
  *  where G is the Givens matrix grot, determined by sigma and gamma. 
  */
  void applyReflector(const givensRot &grot, size_t k, size_t l, EigenVector& p, EigenVector& q);
  

  // @brief Logging device.
  static tarch::logging::Log _log;

  EigenMatrix _Q;
  EigenMatrix _R;
  
  size_t _cols;
  size_t _rows;
  double _omega;
  double _theta;
  double _sigma;
  
};

}}} // namespace precice, cplscheme, impl

#endif /* PRECICE_QRFACTORIZATION_HPP_ */
