
#include "QRFactorization.hpp"
#include "tarch/la/MatrixVectorOperations.h"
#include "utils/Dimensions.hpp"
#include "tarch/la/Scalar.h"

#include <time.h>
#include <math.h> 

namespace precice {
namespace cplscheme {
namespace impl {

tarch::logging::Log QRFactorization::
      _log("precice::cplscheme::impl::QRFactorization");

      

      
/**
 * Constructor
 */
QRFactorization::QRFactorization(
  EigenMatrix Q, 
  EigenMatrix R, 
  size_t rows, 
  size_t cols, 
  double omega, 
  double theta, 
  double sigma)
  :
  _Q(Q),
  _R(R),
  _rows(rows),
  _cols(cols),
  _omega(omega),
  _theta(theta),
  _sigma(sigma)
{
  assertion2(_R.rows() == _cols, _R.rows(), _cols);
  assertion2(_R.cols() == _cols, _R.cols(), _cols);
  assertion2(_Q.cols() == _cols, _Q.cols(), _cols);
  assertion2(_Q.rows() == _rows, _Q.rows(), _rows);
}


/**
 * Constructor
 */
QRFactorization::QRFactorization(
  DataMatrix A, 
  double omega, 
  double theta, 
  double sigma)
  :
  _Q(),
  _R(),
  _rows(A.rows()),
  _cols(0),
  _omega(omega),
  _theta(theta),
  _sigma(sigma)
{
  size_t m = A.cols();
  for (size_t k=0; k<m; k++)
  {
     EigenVector v(m);
     for(size_t i=0; i<_rows; i++)
       v(i) = A(i,k);
     insertColumn(k,v);
  }
  assertion2(_R.rows() == _cols, _R.rows(), _cols);
  assertion2(_R.cols() == _cols, _R.cols(), _cols);
  assertion2(_Q.cols() == _cols, _Q.cols(), _cols);
  assertion2(_Q.rows() == _rows, _Q.rows(), _rows);
  assertion2(_cols == m, _cols, m);
}

/**
 * Constructor
 */
QRFactorization::QRFactorization(
  EigenMatrix A, 
  double omega, 
  double theta, 
  double sigma)
  :
  _Q(),
  _R(),
  _rows(A.rows()),
  _cols(0),
  _omega(omega),
  _theta(theta),
  _sigma(sigma)
{
  size_t m = A.cols();
  for (size_t k=0; k<m; k++)
  {
     EigenVector v = A.col(k);
     //for(int i=0; i<_rows; i++)
     //  v(i) = A(i,k);
     insertColumn(k,v);
  }
  assertion2(_R.rows() == _cols, _R.rows(), _cols);
  assertion2(_R.cols() == _cols, _R.cols(), _cols);
  assertion2(_Q.cols() == _cols, _Q.cols(), _cols);
  assertion2(_Q.rows() == _rows, _Q.rows(), _rows);
  assertion2(_cols == m, _cols, m);
}

      
 
    
/**
 * updates the factorization A=Q[1:n,1:m]R[1:m,1:n] when the kth column of A is deleted. 
 * Returns the deleted column v(1:n)
 */
void QRFactorization::deleteColumn(size_t k)
{
  assertion1(k >= 0, k);
  assertion1(k <= _cols, k);
  
  // maintain decomposition and orthogonalization through application of givens rotations
  double x,y;
  for(size_t l=k; l<_cols-1; l++)
  {
    QRFactorization::givensRot grot;
    computeReflector(grot, _R(l,l+1), _R(l+1,l+1));
    EigenVector Rr1 = _R.row(l);
    EigenVector Rr2 = _R.row(l+1);
    applyReflector(grot, l+2, _cols, Rr1, Rr2);
    _R.row(l) = Rr1;
    _R.row(l+1) = Rr2;
    EigenVector Qc1 = _Q.col(l);
    EigenVector Qc2 = _Q.col(l+1);
    applyReflector(grot, 1, _rows, Qc1, Qc2);
    _Q.col(l) = Qc1;
    _Q.col(l+1) = Qc2;
  }
   
  // copy values and resize R and Q
  for(size_t j=k; k<_cols-1; k++)
  {
    for(size_t i=0; i<=j; i++)
    {
      _R(i,j) = _R(i,j+1);
    }
  }
  _R.conservativeResize(_cols-1, _cols-1);
  //_Q.conservativeResize(Eigen::NoChange_t, _cols-1);
  _Q.conservativeResize(_rows, _cols-1);
  _cols--;
  
  assertion2(_Q.cols() == _cols, _Q.cols(), _cols);
  assertion2(_Q.rows() == _rows, _Q.rows(), _rows);
  assertion2(_R.cols() == _cols, _R.cols(), _cols);
  assertion2(_R.rows() == _cols, _Q.rows(), _cols);
}

      
void QRFactorization::insertColumn(size_t k, DataValues& v)
{
   EigenVector _v(v.size());
   for(size_t i=0; i<v.size();i++)
   {
     _v(i) = v(i);
   }
   insertColumn(k, _v);
}
      
      
void QRFactorization::insertColumn(size_t k, EigenVector& v)
{
  
  assertion1(k >= 0, k);
  assertion1(k <= _cols, k);
  assertion2(v.size() == _rows, v.size(), _rows);
  
  // resize R(1:m, 1:m) -> R(1:m+1, 1:m+1)
  _R.conservativeResize(_cols+1,_cols+1);
  _cols++;
  for(size_t j=_cols-2; j >= k; j--)
  {
   for(size_t i=0; i<j; i++)
   {
     _R(i,j+1) = _R(i,j);
   }
  }
  for(size_t j=k+1; k<_cols; j++)
  {
    _R(j,j) = 0.;
  }
  
  assertion2(_R.cols() == _cols, _R.cols(), _cols);
  assertion2(_R.rows() == _cols, _Q.rows(), _cols);
  
  // orthogonalize v to columns of Q
  EigenVector u(_cols);
  double rho = 0;
  int err = orthogonalize(v, u, rho, _cols-1);
  //_Q.conservativeResize(Eigen::NoChange_t, _cols);
  _Q.conservativeResize(_rows, _cols);
  _Q.col(_cols-1) = v;
  
  assertion2(u(_cols-1) == rho, u.tail(1), rho);
  assertion2(_Q.cols() == _cols, _Q.cols(), _cols);
  assertion2(_Q.rows() == _rows, _Q.rows(), _rows);
  
  // maintain decomposition and orthogonalization through application of givens rotations
  double x,y;
  for(size_t l=_cols-2; l>=k; l--)
  {
    QRFactorization::givensRot grot;
    computeReflector(grot, u(l), u(l+1));
    EigenVector Rr1 = _R.row(l);
    EigenVector Rr2 = _R.row(l+1);
    applyReflector(grot, l+1, _cols, Rr1, Rr2);
    _R.row(l) = Rr1;
    _R.row(l+1) = Rr2;
    EigenVector Qc1 = _Q.col(l);
    EigenVector Qc2 = _Q.col(l+1);
    applyReflector(grot, 1, _rows, Qc1, Qc2);
    _Q.col(l) = Qc1;
    _Q.col(l+1) = Qc2;
  }
  for(size_t i=0; i<=k; i++)
  {
    _R(i,k) = u(i);
  }
}

      
/**
 * @short assuming Q(1:n,1:m) has nearly orthonormal columns, this procedure
 *   orthogonlizes v(1:n) to the columns of Q, and normalizes the result.
 *   r(1:n) is the array of Fourier coefficients, and rho is the distance
 *   from v to range of Q, r and its corrections are computed in double
 *   precision.
 */
int QRFactorization::orthogonalize(
  EigenVector& v, 
  EigenVector& r, 
  double& rho,
  size_t colNum)
{
   bool restart = false;
   bool null = false;
   bool termination = false;
   double rho0 = 0., rho1 = 0.;
   double t = 0;
   EigenVector u = EigenVector::Zero(_rows);
   EigenVector s = EigenVector::Zero(colNum);
   r = EigenVector::Zero(_cols);
   rho = v.norm();
   rho0 = rho;
   size_t k = 0;
   while (!termination)
   {
     // take a gram-schmidt iteration, ignoring r on later steps if previous v was null
     u = EigenVector::Zero(_rows);
     for(size_t j=0; j<colNum; j++)
     {
       double ss = 0;
       for(size_t i=0; i<_rows; i++)
       {
         ss = ss + _Q(i,j) * v(i);
       }
       t = ss;
       s(j) = t;
       for(size_t i=0; i<_rows; i++)
       {
	 u(i) = u(i) + _Q(i,j) * t;
       }
     }
     if(!null)
     {
       for(size_t j=0; j<colNum; j++)
       {
	 r(j) = r(j) + s(j);
       }
     }
     for(size_t i=0; i<_rows; i++)
     {
       v(i) = v(i) - u(i);
     }
     rho1 = v.norm();
     t = s.norm();
     k++;
     
     // treat the special case m=n
     if(_rows == colNum)
     {
      v = EigenVector::Zero(_rows);
      rho = 0.;
      return 0;
     }
     
     // test for nontermination
     if(rho0 + _omega * t >= _theta *rho1)
     {
       // exit to fail of too many iterations
       if (k >= 4)
       {
         std::cout<<"\ntoo many iterations in orthogonalize, termination failed\n";
         return -1;
       }
       if(!restart && rho1 <= rho*_sigma)
       {
	 restart = true;
	 // find first roee of minimal length of Q
	 u = EigenVector::Zero(_rows);
	 for(size_t j=0; j<colNum; j++)
	 {
	   for(size_t i=0; i<_rows; i++)
	   {
	     u(i) = u(i) + _Q(i,j)*_Q(i,j);
	   }
	 }
	 t = 2;
	 for(size_t i=0; i<_rows; i++)
	 {
	   if(u(i) < t)
	   {
	     k = i;
	     t = u(k);
	   }
	 }
	 
	 // take correct action if v is null
	 if(rho1 == 0)
	 {
	   null = true;
	   rho1 = 1;
	 }
	 // reinitialize v and k
	 v = EigenVector::Zero(_rows);
	 v(k) = rho1;
	 k = 0;
       }
       rho0 = rho1;
     }else
     {
       termination = true;
     }
   }
   
   // normalize v
   v /= rho1;
   rho = null ? 0 : rho1;
   r(colNum) = rho;
   return 0;
}      

   
      

/**
 * @short computes parameters for givens matrix G for which  (x,y)G = (z,0). replaces (x,y) by (z,0)
 */
void QRFactorization::computeReflector(
  QRFactorization::givensRot& grot, 
  double& x, 
  double& y)
{
  double u = x;
  double v = y;
  if (v==0)
  {
    grot.sigma = 0;
    grot.gamma = 1;
  }else
  {
    double mu = std::max(std::fabs(u),std::fabs(v));
    double t = mu * std::sqrt(std::pow(u/mu, 2) + std::pow(v/mu, 2));
    t *= (u < 0) ? -1 : 1;
    grot.gamma = u/t;
    grot.sigma = v/t;
    x = t;
    y = 0;
  }
}

      
/**
 *  @short this procedure replaces the two column matrix [p(k:l-1), q(k:l-1)] by [p(k:l), q(k:l)]*G, 
 *  where G is the Givens matrix grot, determined by sigma and gamma. 
 */
void QRFactorization::applyReflector(
  const QRFactorization::givensRot& grot, 
  size_t k, 
  size_t l, 
  EigenVector& p, 
  EigenVector& q)
{
  double nu = grot.sigma/(1.+grot.gamma);
  for(size_t j=k; j<l; j++)
  {
    double u = p(j);
    double v = q(j);
    double t = u*grot.gamma + v*grot.sigma;
    p(j) = t;
    q(j) = (t + u) * nu - v;
  }
}


QRFactorization::EigenMatrix& QRFactorization::matrixQ()
{
  return _Q;
}

QRFactorization::EigenMatrix& QRFactorization::matrixR()
{
  return _R;
}

size_t QRFactorization::cols()
{
  return _cols;
}

size_t QRFactorization::rows()
{
  return _rows;
}

void QRFactorization::reset()
{
  _Q.resize(0,0);
  _R.resize(0,0);
  _cols = 0;
  _rows = 0;
}

void QRFactorization::reset(
  EigenMatrix Q, 
  EigenMatrix R, 
  size_t rows, 
  size_t cols, 
  double omega, 
  double theta, 
  double sigma)
{
  _Q = Q;
  _R = R;
  _rows = rows;
  _cols = cols;
  _omega = omega;
  _theta = theta;
  _sigma = sigma;
  assertion2(_R.rows() == _cols, _R.rows(), _cols);
  assertion2(_R.cols() == _cols, _R.cols(), _cols);
  assertion2(_Q.cols() == _cols, _Q.cols(), _cols);
  assertion2(_Q.rows() == _rows, _Q.rows(), _rows);
}

void QRFactorization::reset(
  EigenMatrix A, 
  double omega, 
  double theta, 
  double sigma)
{
  _Q.resize(0,0);
  _R.resize(0,0);
  _cols = 0;
  _rows = 0;
  _omega = omega;
  _theta = theta;
  _sigma = sigma;
  
  size_t m = A.cols();
  for (size_t k=0; k<m; k++)
  {
     EigenVector v = A.col(k);
     //for(int i=0; i<_rows; i++)
     //  v(i) = A(i,k);
     insertColumn(k,v);
  }
  assertion2(_R.rows() == _cols, _R.rows(), _cols);
  assertion2(_R.cols() == _cols, _R.cols(), _cols);
  assertion2(_Q.cols() == _cols, _Q.cols(), _cols);
  assertion2(_Q.rows() == _rows, _Q.rows(), _rows);
  assertion2(_cols == m, _cols, m);
}



}}} // namespace precice, cplscheme, impl