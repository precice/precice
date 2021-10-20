#pragma once

#include <Eigen/Core>
#include <fstream>
#include <limits>
#include <string>
#include <vector>
#include "logging/Logger.hpp"
#include "mesh/SharedPointer.hpp"

namespace precice {
namespace acceleration {
namespace impl {

/**
 * @brief Class that provides functionality for a dynamic QR-decomposition, that can be updated
 * in O(mn) flops if a column is inserted or deleted.
 * The new colmn is orthogonalized to the existing columns in Q using a modified GramSchmidt algorithm.
 * The zero-elements are generated using suitable givens-roatations.
 * The Interface provides fnctions such as insertColumn, deleteColumn at arbitrary position an push or pull
 * column at front or back, resp.
 */
class QRFactorization {
public:
  /**
   * @brief Constructor.
   * @param theta - singularity limit for reothogonalization ||v_orth|| / ||v|| <= 1/theta
   */
  QRFactorization(
      int    filter = 0,
      double omega  = 0,
      double theta  = 1. / 0.7,
      double sigma  = std::numeric_limits<double>::min());

  /**
   * @brief Constructor.
   * @param theta - singularity limit for reothogonalization ||v_orth|| / ||v|| <= 1/theta
   */
  QRFactorization(
      Eigen::MatrixXd A,
      int             filter,
      double          omega = 0,
      double          theta = 1. / 0.7,
      double          sigma = std::numeric_limits<double>::min());

  /**
   * @brief Constructor.
   * @param theta - singularity limit for reothogonalization ||v_orth|| / ||v|| <= 1/theta
   */
  QRFactorization(
      Eigen::MatrixXd Q,
      Eigen::MatrixXd R,
      int             rows,
      int             cols,
      int             filter,
      double          omega = 0,
      double          theta = 1. / 0.7,
      double          sigma = std::numeric_limits<double>::min());

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
      Eigen::MatrixXd const &Q,
      Eigen::MatrixXd const &R,
      int                    rows,
      int                    cols,
      double                 omega = 0,
      double                 theta = 1. / 0.7,
      double                 sigma = std::numeric_limits<double>::min());

  /**
    * @brief resets the QR factorization to be the factorization of A = QR
    */
  void reset(
      Eigen::MatrixXd const &A,
      int                    globalRows,
      double                 omega = 0,
      double                 theta = 1. / 0.7,
      double                 sigma = std::numeric_limits<double>::min());

  /**
    * @brief inserts a new column at arbitrary position and updates the QR factorization
    * This function works on the memory of v, thus changes the Vector v.
    */
  bool insertColumn(int k, const Eigen::VectorXd &v, double singularityLimit = 0);

  /**
   * @brief updates the factorization A=Q[1:n,1:m]R[1:m,1:n] when the kth column of A is deleted.
   * Returns the deleted column v(1:n)
   */
  void deleteColumn(int k);

  /**
    * @brief inserts a new column at position 0, i.e., shifts right and inserts at first position
    * and updates the QR factorization.
    * This function works on the memory of v, thus changes the Vector v.
    */
  void pushFront(const Eigen::VectorXd &v);

  /**
    * @brief inserts a new column at position _cols-1, i.e., appends a column at the end
    * and updates the QR factorization
    * This function works on the memory of v, thus changes the Vector v.
    */
  void pushBack(const Eigen::VectorXd &v);

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
    * @brief filters the least squares system, i.e., the decomposition Q*R = V according
    * to the defined filter technique. This is done to ensure good conditioning
    * @param [out] delIndices - a vector of indices of deleted columns from the LS-system
    */
  void applyFilter(double singularityLimit, std::vector<int> &delIndices, Eigen::MatrixXd &V);

  /**
    * @brief returns a matrix representation of the orthogonal matrix Q
    */
  Eigen::MatrixXd &matrixQ();

  /**
    * @brief returns a matrix representation of the upper triangular matrix R
    */
  Eigen::MatrixXd &matrixR();

  // @brief returns the number of columns in the QR-decomposition
  int cols() const;
  // @brief returns the number of rows in the QR-decomposition
  int rows() const;

  // @brief optional file-stream for logging output
  void setfstream(std::fstream *stream);

  // @brief set number of global rows for the master-slave case
  void setGlobalRows(int gr);

  // @brief sets the filtering technique to maintain good conditioning of the least squares system
  void setFilter(int filter);

private:
  struct givensRot {
    int    i, j;
    double sigma, gamma;
  };

  /**
  * @short assuming Q(1:n,1:m) has nearly orthonormal columns, this procedure
  *   orthogonlizes v(1:n) to the columns of Q, and normalizes the result.
  *   r(1:n) is the array of Fourier coefficients, and rho is the distance
  *   from v to range of Q, r and its corrections are computed in double
  *   precision.
  *   This method tries to re-orthogonalize the matrix Q for a maximum of 4 iterations
  *   if ||v_orth|| / ||v|| <= 1/theta is toot small, i.e. the gram schmidt process is iterated.
  *   If ||v_orth|| / ||v|| <= std::numeric_limit, a unit vector that is orthogonal to Q is inserted
  *   and rho is set to 0. i.e., R has a zero on the diagonal in the respective column.
  */
  int orthogonalize_stable(Eigen::VectorXd &v, Eigen::VectorXd &r, double &rho, int colNum);

  /**
   * @short assuming Q(1:n,1:m) has nearly orthonormal columns, this procedure
  *   orthogonlizes v(1:n) to the columns of Q, and normalizes the result.
  *   r(1:n) is the array of Fourier coefficients, and rho is the distance
  *   from v to range of Q, r and its corrections are computed in double
  *   precision.
  *   This method tries to re-orthogonalize the matrix Q for a maximum of 4 iterations
  *   if ||v_orth|| / ||v|| <= 1/theta is toot small, i.e. the gram schmidt process is iterated.
  *
  *   Difference to the method orthogonalize_stable():
  *   if ||v_orth||/||v|| approx 0, no unit vector is inserted.
   */
  int orthogonalize(Eigen::VectorXd &v, Eigen::VectorXd &r, double &rho, int colNum);

  /**
  * @short computes parameters for givens matrix G for which  (x,y)G = (z,0). replaces (x,y) by (z,0)
  */
  void computeReflector(givensRot &grot, double &x, double &y);

  /**
  *  @short this procedure replaces the two column matrix [p(k:l-1), q(k:l-1)] by [p(k:l), q(k:l)]*G,
  *  where G is the Givens matrix grot, determined by sigma and gamma.
  */
  void applyReflector(const givensRot &grot, int k, int l, Eigen::VectorXd &p, Eigen::VectorXd &q);

  logging::Logger _log{"acceleration::QRFactorization"};

  Eigen::MatrixXd _Q;
  Eigen::MatrixXd _R;

  int _rows;
  int _cols;

  int    _filter;
  double _omega;
  double _theta;
  double _sigma;

  // @brief optional infostream that writes information to file
  std::fstream *_infostream;
  bool          _fstream_set;

  int _globalRows;
};

} // namespace impl
} // namespace acceleration
} // namespace precice
