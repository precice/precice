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
 * @brief Class that provides functionality for a dynamic QR-decomposition, that can be updated 
 * in O(mn) flops if a column is inserted or deleted. 
 * The new colmn is orthogonalized to the existing columns in Q using a modified GramSchmidt algorithm.
 * The zero-elements are generated using suitable givens-roatations.
 * The Interface provides fnctions such as insertColumn, deleteColumn at arbitrary position an push or pull 
 * column at front or back, resp. 
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
      double omega=0,
      double theta=1./0.7,
      double sigma=std::numeric_limits<double>::min()
 		    );
  
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
      int rows,
      int cols,
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
	int rows, 
	int cols, 
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
   
   /**
    * @brief resets the QR factorization to be the factorization of A = QR
    */
   void reset(
	DataMatrix A,
	double omega=0,
        double theta=1./0.7,
        double sigma=std::numeric_limits<double>::min());
   
   /**
    * @brief inserts a new column at arbitrary position and updates the QR factorization
    */
   void insertColumn(int k, EigenVector& v);
   void insertColumn(int k, DataValues& v);
   
   /**
   * @brief updates the factorization A=Q[1:n,1:m]R[1:m,1:n] when the kth column of A is deleted. 
   * Returns the deleted column v(1:n)
   */
   void deleteColumn(int k);
   
   /**
    * @brief inserts a new column at position 0, i.e., shifts right and inserts at first position
    * and updates the QR factorization
    */
   void pushFront(EigenVector& v);
   
   /**
    * @brief inserts a new column at position _cols-1, i.e., appends a column at the end
    * and updates the QR factorization
    */
   void pushBack(EigenVector& v);
   
   /**
    * @brief inserts a new column at position 0, i.e., shifts right and inserts at first position
    * and updates the QR factorization
    */
   void pushFront(DataValues& v);
   
   /**
    * @brief inserts a new column at position _cols-1, i.e., appends a column at the end
    * and updates the QR factorization
    */
   void pushBack(DataValues& v);
   
   /**
    * @brief deletes the column at position 0, i.e., deletes and shifts columns to the left
    * and updates the QR factorization
    */
   void popFront();
   
   /**
    * @brief deletes the column at position _cols-1, i.e., deletes the last column
    * and updates the QR factorization
    */
   void popBack();
   
   /**
    * @brief returns a matrix representation of the orthogonal matrix Q
    */
   EigenMatrix& matrixQ();

   /**
    * @brief returns a matrix representation of the upper triangular matrix R
    */
   EigenMatrix& matrixR();
   
   // @brief returns the number of columns in the QR-decomposition
   int cols();
   // @brief returns the number of rows in the QR-decomposition
   int rows();

   // @brief optional file-stream for logging output
   void setfstream(std::fstream* stream);

   // @brief set number of global rows for the master-slave case
   void setGlobalRows(int gr);

private:
  
  struct givensRot{
    int i, j;
    double sigma, gamma;
  };
   
  /**
  * @short assuming Q(1:n,1:m) has nearly orthonormal columns, this procedure
  *   orthogonlizes v(1:n) to the columns of Q, and normalizes the result.
  *   r(1:n) is the array of Fourier coefficients, and rho is the distance
  *   from v to range of Q, r and its corrections are computed in double
  *   precision.
  */
  int orthogonalize(EigenVector& v, EigenVector& r, double &rho, int colNum);
  
  /**
  * @short computes parameters for givens matrix G for which  (x,y)G = (z,0). replaces (x,y) by (z,0)
  */
  void computeReflector(givensRot &grot, double &x, double &y);
  
  /**
  *  @short this procedure replaces the two column matrix [p(k:l-1), q(k:l-1)] by [p(k:l), q(k:l)]*G, 
  *  where G is the Givens matrix grot, determined by sigma and gamma. 
  */
  void applyReflector(const givensRot &grot, int k, int l, EigenVector& p, EigenVector& q);
  

  // @brief Logging device.
  static tarch::logging::Log _log;

  EigenMatrix _Q;
  EigenMatrix _R;
  
  int _rows;
  int _cols;

  double _omega;
  double _theta;
  double _sigma;
  
  // @brief optional infostream that writes information to file
  std::fstream* _infostream;
  bool _fstream_set;

  int _globalRows;

};

}}} // namespace precice, cplscheme, impl

#endif /* PRECICE_QRFACTORIZATION_HPP_ */
