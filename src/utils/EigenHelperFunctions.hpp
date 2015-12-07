/*
 * EigenHelperFunctions.hpp
 *
 *  Created on: Dec 7, 2015
 *      Author: scheufks
 */

#ifndef EIGENHELPERFUNCTIONS_HPP_
#define EIGENHELPERFUNCTIONS_HPP_

#include "Eigen/Dense"

namespace precice {
namespace utils {


void shiftSetFirst
(
    Eigen::MatrixXd& A, Eigen::VectorXd& v)
{
  assertion2(v.size() == A.rows(), v.size(), A.rows());
  int n = A.rows(), m = A.cols();
  //A.bottomRightCorner(n, m - 1) = A.topLeftCorner(n, m - 1);
  for(auto i = A.cols()-1; i > 0; i--)
        A.col(i) = A.col(i-1);
  A.col(0) = v;
}

void appendFront
(
    Eigen::MatrixXd& A, Eigen::VectorXd& v)
{
  int n = A.rows(), m = A.cols();
  if (n <= 0 && m <= 0) {
    A = v;
  } else {
    assertion2(v.size() == n, v.size(), A.rows());
    A.conservativeResize(n, m + 1);
    //A.topRightCorner(n, m) = A.topLeftCorner(n, m); // bad error, reason unknown!
    for(auto i = A.cols()-1; i > 0; i--)
      A.col(i) = A.col(i-1);
    A.col(0) = v;
  }
}

void append
(
    Eigen::MatrixXd& A, Eigen::MatrixXd& B)
{
  int n = A.rows(), m = A.cols();
  if (n <= 0 && m <= 0) {
    A = B;
  } else {
    assertion2(B.rows() == n, B.rows(), A.rows());
    A.conservativeResize(n, m + B.cols());
    for(auto i = m; i < m + B.cols(); i++)
      A.col(i) = B.col(i-m);
  }
}

void append
(
    Eigen::MatrixXd& A, Eigen::VectorXd& v)
{
  int n = A.rows(), m = A.cols();
  if (n <= 0 && m <= 0) {
    A = v;
  } else {
    assertion2(v.size() == n, v.size(), A.rows());
    A.conservativeResize(n, m + 1);
    A.col(m+1) = v;
  }
}

void append
(
    Eigen::VectorXd& v, Eigen::VectorXd& app)
{
  int n = v.size();
  v.conservativeResize(n + app.size());
  for(int i = 0; i < app.size(); i++)
    v(n+i) = app(i);
}

void removeColumnFromMatrix
(
    Eigen::MatrixXd& A, int col)
{
  assertion2(col < A.cols() && col >= 0, col, A.cols())
  for (int j = col; j < A.cols() - 1; j++)
    A.col(j) = A.col(j + 1);

  A.conservativeResize(A.rows(), A.cols() - 1);
}



}}
#endif /* EIGENHELPERFUNCTIONS_HPP_ */
