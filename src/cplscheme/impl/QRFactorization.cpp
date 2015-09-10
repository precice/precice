
#include "QRFactorization.hpp"
#include "tarch/la/MatrixVectorOperations.h"
#include "utils/Dimensions.hpp"
#include "utils/Globals.hpp"
#include "tarch/la/Scalar.h"
#include "cplscheme/impl/BaseQNPostProcessing.hpp"

#include <iostream>
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
  int rows, 
  int cols, 
  int filter,
  double omega, 
  double theta, 
  double sigma)

  :
  _Q(Q),
  _R(R),
  _rows(rows),
  _cols(cols),
  _filter(filter),
  _omega(omega),
  _theta(theta),
  _sigma(sigma),
  _infostream(),
  _fstream_set(false),
  _globalRows(rows)
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
  int filter,
  double omega, 
  double theta, 
  double sigma)
  :
  _Q(),
  _R(),
  _rows(A.rows()),
  _cols(0),
  _filter(filter),
  _omega(omega),
  _theta(theta),
  _sigma(sigma),
  _infostream(),
  _fstream_set(false),
  _globalRows(A.rows())
{
  int m = A.cols();
  for (int k=0; k<m; k++)
  {
     EigenVector v(_rows);
     for(int i=0; i<_rows; i++)
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
  int filter,
  double omega, 
  double theta, 
  double sigma)
  :
  _Q(),
  _R(),
  _rows(A.rows()),
  _cols(0),
  _filter(filter),
  _omega(omega),
  _theta(theta),
  _sigma(sigma),
  _infostream(),
  _fstream_set(false),
  _globalRows(A.rows())
{
  int m = A.cols();
  for (int k=0; k<m; k++)
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
 * Constructor
 */
QRFactorization::QRFactorization(
  int filter,
  double omega, 
  double theta, 
  double sigma)
  :
  _Q(),
  _R(),
  _rows(0),
  _cols(0),
  _filter(filter),
  _omega(omega),
  _theta(theta),
  _sigma(sigma),
  _infostream(),
  _fstream_set(false),
  _globalRows(0)
{}

void QRFactorization::applyFilter(double singularityLimit, std::vector<int>& delIndices, DataMatrix& V)
{
	EigenMatrix _V(V.rows(), V.cols());
	for(int i = 0; i < _V.rows(); i++)
		for(int j = 0; j < _V.cols(); j++){
			_V(i,j) = V(i,j);
		}
	applyFilter(singularityLimit, delIndices, _V);
}
      
void QRFactorization::applyFilter(double singularityLimit, std::vector<int>& delIndices, EigenMatrix& V)
{
	preciceTrace("applyFilter()");
	delIndices.resize(0);
	if(_filter == BaseQNPostProcessing::QR1FILTER || _filter == BaseQNPostProcessing::QR1FILTER_ABS)
	{
		bool linearDependence = true;
		std::vector<int> delFlag(_cols, 0);
		int delCols = 0;
		while (linearDependence) {
			linearDependence = false;
			int index = 0; // actual index of checked column, \in [0, _cols] and _cols is decreasing
			if(_cols > 1){
				for (int i = 0; i < delFlag.size(); i++) {
					// index is not incremented, if columns has been deleted in previous rounds
					if(delFlag[i] > 0) continue;

					// QR1-filter
					if(index >= cols()) break;
					assertion2(index < _cols, index, _cols);
					double factor = (_filter == BaseQNPostProcessing::QR1FILTER_ABS) ? 1.0 : _R.norm();
					if (std::fabs(_R(index, index)) < singularityLimit * factor) {

						linearDependence = true;
						deleteColumn(index);
						delFlag[i]++;
						delIndices.push_back(i);
						delCols++;
						//break;
						index--;  	// check same column index, as cols are shifted left
					}
					assertion2(delCols+_cols == delFlag.size(), (delCols+_cols), delFlag.size());
					index++;
				}
			}
		}
	}else if(_filter == BaseQNPostProcessing::QR2FILTER)
	{
		  _Q.resize(0,0);
		  _R.resize(0,0);
		  _cols = 0;
		  _rows = V.rows();
		  // starting with the most recent input/output information, i.e., the latest column
		  // which is at position 0 in _matrixV (latest information is never filtered out!)
		  for (int k=0; k < V.cols(); k++)
		  {
		     EigenVector v = V.col(k);
		     bool inserted = insertColumn(k, v, singularityLimit);
		     if (!inserted){
		    	 delIndices.push_back(k);
		     }
		  }
	}
}
 
    
/**
 * updates the factorization A=Q[1:n,1:m]R[1:m,1:n] when the kth column of A is deleted. 
 * Returns the deleted column v(1:n)
 */
void QRFactorization::deleteColumn(int k)
{

  preciceTrace("deleteColumn()");

  assertion1(k >= 0, k);
  assertion2(k < _cols, k, _cols);
  
  // maintain decomposition and orthogonalization by application of givens rotations
  
  for(int l=k; l<_cols-1; l++)
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
    applyReflector(grot, 0, _rows, Qc1, Qc2);
    _Q.col(l) = Qc1;
    _Q.col(l+1) = Qc2;
  }
  // copy values and resize R and Q
  for(int j=k; j<_cols-1; j++)
  {
    for(int i=0; i<=j; i++)
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

      
bool QRFactorization::insertColumn(int k, DataValues& v, double singularityLimit)
{
   EigenVector _v(v.size());
   for(int i=0; i<v.size();i++)
   {
     _v(i) = v(i);
   }
   return insertColumn(k, _v, singularityLimit);
}
      
      
bool QRFactorization::insertColumn(int k, EigenVector& v, double singularityLimit)
{
  preciceTrace("insertColumn()");

  if(_cols == 0)
    _rows = v.size();
  
  bool applyFilter = (singularityLimit > 0.0);

  assertion1(k >= 0, k);
  assertion1(k <= _cols, k);
  assertion2(v.size() == _rows, v.size(), _rows);
  
  _cols++;

  // orthogonalize v to columns of Q
  EigenVector u(_cols);
  double rho_orth = 0., rho0 = 0.;
  if(applyFilter) rho0 = utils::MasterSlave::l2norm(v);

  int err = orthogonalize(v, u, rho_orth, _cols-1);
  // if err < 0, the column cannot be orthogonalized to the system
  // if rho_orth = 0, either v = 0 or the system is quadratic.
  // Discard column v in all cases.
  if(err < 0 || rho_orth == 0.0){
	  preciceDebug("discarding column because err < 0 or rho_orth == 0.0, i.e., too many iterations in orthogonalize or quadratic system. err: "<<err<<", rho_orth: "<<rho_orth);
	  _cols--;
	  return false;
  }

  // QR2-filter based on the new information added to the orthogonal system.
  // if the new column incorporates less new information to the system than a
  // prescribed threshold, the column is discarded
  // rho_orth: the norm of the orthogonalized (but not normalized) column
  // rho0:     the norm of the initial column that is to be inserted
  if(applyFilter && (rho0 * singularityLimit  > rho_orth)){
	preciceDebug("discarding column as it is filtered out by the QR2-filter: rho0*eps > rho_orth: "<<rho0*singularityLimit<<" > "<<rho_orth);
    _cols--;
    return false;
  }


  // resize R(1:m, 1:m) -> R(1:m+1, 1:m+1)
  _R.conservativeResize(_cols,_cols);
  _R.col(_cols-1) = EigenVector::Zero(_cols);
  _R.row(_cols-1) = EigenVector::Zero(_cols);
  
  for(int j=_cols-2; j >= k; j--)
  {
    for(int i=0; i<=j; i++)
    {
      _R(i,j+1) = _R(i,j);
    }
  }
  
  for(int j=k+1; j<_cols; j++)
  {
    _R(j,j) = 0.;
  }
  
  assertion2(_R.cols() == _cols, _R.cols(), _cols);
  assertion2(_R.rows() == _cols, _R.rows(), _cols);
  
  // resize Q(1:n, 1:m) -> Q(1:n, 1:m+1)
  _Q.conservativeResize(_rows, _cols);
  _Q.col(_cols-1) = v;
  
  assertion2(_Q.cols() == _cols, _Q.cols(), _cols);
  assertion2(_Q.rows() == _rows, _Q.rows(), _rows);
  
  // maintain decomposition and orthogonalization by application of givens rotations
  for(int l=_cols-2; l>=k; l--)
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
    applyReflector(grot, 0, _rows, Qc1, Qc2);
    _Q.col(l) = Qc1;
    _Q.col(l+1) = Qc2;
  }
  for(int i=0; i<=k; i++)
  {
    _R(i,k) = u(i);
  }

  return true;
}

      
/**
 * @short assuming Q(1:n,1:m) has nearly orthonormal columns, this procedure
 *   orthogonlizes v(1:n) to the columns of Q, and normalizes the result.
 *   r(1:n) is the array of Fourier coefficients, and rho is the distance
 *   from v to range of Q, r and its corrections are computed in double
 *   precision.
 *
 *   @return Returns the number of gram-schmidt iterations needed to orthogobalize the
 *   new vector to the existing system. If more then 4 iterations were needed, -1 is
 *   returned and the new column should not be inserted into the system.
 */
int QRFactorization::orthogonalize(
  EigenVector& v, 
  EigenVector& r, 
  double& rho,
  int colNum)
{
   preciceTrace("orthogonalize()");

   // serial case
   if(not utils::MasterSlave::_masterMode && not utils::MasterSlave::_slaveMode)
   {
	   assertion2(_globalRows == _rows, _globalRows, _rows);
   // master-slave case
   }else{
	   assertion3(_globalRows != _rows, _globalRows, _rows, utils::MasterSlave::_rank);
   }

   bool restart = false;
   bool null = false;
   bool termination = false;
   double rho0 = 0., rho1 = 0.;
   double t = 0;
   EigenVector u = EigenVector::Zero(_rows);
   EigenVector s = EigenVector::Zero(colNum);
   r = EigenVector::Zero(_cols);
   
   rho = utils::MasterSlave::l2norm(v); // distributed l2norm
   rho0 = rho;
   int k = 0;
	while (!termination) {
		// take a gram-schmidt iteration, ignoring r on later steps if previous v was null
		u = EigenVector::Zero(_rows);
		for (int j = 0; j < colNum; j++) {

			/*
			 * dot-product <_Q(:,j), v >
			 */
			EigenVector Qc = _Q.col(j);

			// dot product <_Q(:,j), v> =: r_ij
			double ss = utils::MasterSlave::dot(Qc, v);
			t = ss;
			// save r_ij in s(j) = column of R
			s(j) = t;
			// u is the sum of projections r_ij * _Q(i,:) =  _Q(i,:) * <_Q(:,j), v>
			for (int i = 0; i < _rows; i++) {
				u(i) = u(i) + _Q(i, j) * t;
			}
		}
		if (!null) {
			// add over all runs: r_ij = r_ij_prev + r_ij
			for (int j = 0; j < colNum; j++) {
				r(j) = r(j) + s(j);
			}
		}
		// subtract projections from v, v is now orthogonal to columns of _Q
		for (int i = 0; i < _rows; i++) {
			v(i) = v(i) - u(i);
		}
		// rho1 = norm of orthogonalized new column v_tilde (though not normalized)
		rho1 = utils::MasterSlave::l2norm(v); // distributed l2norm

		// t = norm of r_(:,j) with j = colNum-1
		t = utils::MasterSlave::l2norm(s); // distributed l2norm
		k++;

		// treat the special case m=n
		// Attention (Master-Slave): Here, we need to compare the global _rows with colNum and NOT the local
		// rows on the processor.
		if (_globalRows == colNum) {
			preciceWarning("orthogonalize()", "The least-squares system matrix is quadratic, i.e., the new column cannot be orthogonalized (and thus inserted) to the LS-system.\nOld columns need to be removed.");
			v = EigenVector::Zero(_rows);
			rho = 0.;
			return k;
		}

		/**   - test for nontermination -
		 *  rho0 = |v_init|, t = |r_(i,cols-1)|, rho1 = |v_orth|
		 *  rho1 is small, if the new information incorporated in v is small,
		 *  i.e., the part of v orthogonal to _Q is small.
		 *  if rho1 is very small it is possible, that we are adding (more or less)
		 *  only round-off errors to the decomposition. Later normalization will scale
		 *  this new information so that it is equally weighted as the columns in Q.
		 *  To keep a good orthogonality, some effort is done if comparatively little
		 *  new information is added.
		 */
		if (rho0 + _omega * t >= _theta * rho1) {
			// exit to fail of too many iterations
			if (k >= 4) {
				std::cout
						<< "\ntoo many iterations in orthogonalize, termination failed\n";
				preciceWarning("orthogonalize()", "Matrix Q is not sufficiently orthogonal. Failed to rorthogonalize new column after 4 iterations. New column will be discarded. The least-squares system is very bad conditioned and the quasi-Newton will most probably fail to converge.");
				return -1;
			}
			if (!restart && rho1 <= rho * _sigma) {
			  preciceWarning("orthogonalize()", "Matrix Q is not sufficiently orthogonal, most probably due to bad conditioning of the least squares system. Try to re-orthogonalize the QR system");
				//preciceDebug("[QR-dec] - reorthogonalization");
				if (_fstream_set)
					(*_infostream) << "[QR-dec] - reorthogonalization"
							<< std::endl;

				restart = true;

				/**  - find first row of minimal length of Q -
				 *  the squared l2-norm of each row is computed. Find row of minimal length.
				 *  Start with a new vector v that is zero except for v(k) = rho1, where
				 *  k is the index of the row of Q with minimal length.
				 *  Note: the new information from v is discarded. Q is made orthogonal
				 *        as good as possible.
				 */
				u = EigenVector::Zero(_rows);
				for (int j = 0; j < colNum; j++) {
					for (int i = 0; i < _rows; i++) {
						u(i) = u(i) + _Q(i, j) * _Q(i, j);
					}
				}
				t = 2;

				// ATTENTION: maybe in the following is something wrong in master-slave mode

				for (int i = 0; i < _rows; i++) {
					if (u(i) < t) {
						k = i;
						t = u(k);
					}
				}

				int global_k = k;
				int local_k = 0;
				double local_uk = 0.;
				double global_uk = 0.;
				int rank = 0;

				if (utils::MasterSlave::_slaveMode) {
					utils::MasterSlave::_communication->send(k, 0);
					utils::MasterSlave::_communication->send(u(k), 0);
				}

				if (utils::MasterSlave::_masterMode) {
					global_uk = u(k);
					for (int rankSlave = 1; rankSlave < utils::MasterSlave::_size; rankSlave++) {
						utils::MasterSlave::_communication->receive(local_k, rankSlave);
						utils::MasterSlave::_communication->receive(local_uk, rankSlave);
						if (local_uk < global_uk) {
							rank = rankSlave;
							global_uk = local_uk;
							global_k = local_k;
						}
					}
					if (_fstream_set)
						(*_infostream) << "           global u(k):" << global_uk
								<< ",  global k: " << global_k << ",  rank: "
								<< rank << std::endl;
				}

				// take correct action if v is null
				if (rho1 == 0) {
					null = true;
					rho1 = 1;
				}
				// reinitialize v and k
				v = EigenVector::Zero(_rows);

				// insert rho1 at position k with smallest u(i) = Q(i,:) * Q(i,:)
				if (not utils::MasterSlave::_masterMode && not utils::MasterSlave::_slaveMode) {
					v(k) = rho1;
				} else {
					if (utils::MasterSlave::_rank == rank)
						v(global_k) = rho1;
				}
				k = 0;
			}
			rho0 = rho1;

			// termination, i.e., (rho0 + _omega * t < _theta *rho1)
		} else {
			termination = true;
		}
	}
   
   // normalize v
   v /= rho1;
   rho = null ? 0 : rho1;
   r(colNum) = rho;
   return k;
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
  int k, 
  int l, 
  EigenVector& p, 
  EigenVector& q)
{
  double nu = grot.sigma/(1.+grot.gamma);
  for(int j=k; j<l; j++)
  {
    double u = p(j);
    double v = q(j);
    double t = u*grot.gamma + v*grot.sigma;
    p(j) = t;
    q(j) = (t + u) * nu - v;
  }
}


void QRFactorization::setGlobalRows( int gr)
{
  _globalRows = gr;
}


QRFactorization::EigenMatrix& QRFactorization::matrixQ()
{
  return _Q;
}

QRFactorization::EigenMatrix& QRFactorization::matrixR()
{
  return _R;
}

int QRFactorization::cols()
{
  return _cols;
}

int QRFactorization::rows()
{
  return _rows;
}

void QRFactorization::reset()
{
  _Q.resize(0,0);
  _R.resize(0,0);
  _cols = 0;
  _rows = 0;
  _globalRows = 0;
}

void QRFactorization::reset(
  EigenMatrix Q, 
  EigenMatrix R, 
  int rows, 
  int cols, 
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
  _globalRows = _rows;
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
  _rows = A.rows();
  _omega = omega;
  _theta = theta;
  _sigma = sigma;
  _globalRows = _rows;
  
  int m = A.cols();
  for (int k=0; k<m; k++)
  {
     EigenVector v = A.col(k);
     insertColumn(k,v);
  }
  assertion2(_R.rows() == _cols, _R.rows(), _cols);
  assertion2(_R.cols() == _cols, _R.cols(), _cols);
  assertion2(_Q.cols() == _cols, _Q.cols(), _cols);
  assertion2(_Q.rows() == _rows, _Q.rows(), _rows);
  assertion2(_cols == m, _cols, m);
}

void QRFactorization::reset(
  DataMatrix A, 
  double omega, 
  double theta, 
  double sigma)
{
  _Q.resize(0,0);
  _R.resize(0,0);
  _cols = 0;
  _rows = A.rows();
  _omega = omega;
  _theta = theta;
  _sigma = sigma;
  _globalRows = _rows;
 
  int m = A.cols();
  for (int k=0; k<m; k++)
  {
     EigenVector v(_rows);
     for(int i=0; i<_rows; i++)
       v(i) = A(i,k);
     insertColumn(k,v);
  }
  assertion2(_R.rows() == _cols, _R.rows(), _cols);
  assertion2(_R.cols() == _cols, _R.cols(), _cols);
  assertion2(_Q.cols() == _cols, _Q.cols(), _cols);
  assertion2(_Q.rows() == _rows, _Q.rows(), _rows);
  assertion2(_cols == m, _cols, m);
}

void QRFactorization::pushFront(EigenVector& v)
{
  insertColumn(0, v);
}

void QRFactorization::pushBack(EigenVector& v)
{
  insertColumn(_cols, v);
}

void QRFactorization::pushFront(DataValues& v)
{
  EigenVector _v(v.size());
  for(int i = 0; i<v.size(); i++)
    _v(i) = v(i);
  insertColumn(0, _v);
}

void QRFactorization::pushBack(DataValues& v)
{
  EigenVector _v(v.size());
  for(int i = 0; i<v.size(); i++)
    _v(i) = v(i);
  insertColumn(_cols, _v);
}

void QRFactorization::popFront()
{
  deleteColumn(0);
}

void QRFactorization::popBack()
{
  deleteColumn(_cols-1);
}

void QRFactorization::setfstream(std::fstream* stream)
{
  _infostream = stream;
  _fstream_set = true;
}

void QRFactorization::setFilter(int filter){
	_filter = filter;
}







}}} // namespace precice, cplscheme, impl
