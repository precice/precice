/*
 * EigenHelperFunctions.cpp
 *
 *  Created on: Dec 9, 2015
 *      Author: scheufks
 */

#include "EigenHelperFunctions.hpp"

namespace precice::utils {

void shiftSetFirst(
    Eigen::MatrixXd &A, const Eigen::VectorXd &v)
{
  PRECICE_ASSERT(v.size() == A.rows(), v.size(), A.rows());
  // A.bottomRightCorner(n, m - 1) = A.topLeftCorner(n, m - 1);
  for (auto i = A.cols() - 1; i > 0; i--)
    A.col(i) = A.col(i - 1);
  A.col(0) = v;
}

void appendFront(
    Eigen::MatrixXd &A, Eigen::VectorXd &v)
{
  int n = A.rows(), m = A.cols();
  if (n <= 0 && m <= 0) {
    A = v;
  } else {
    PRECICE_ASSERT(v.size() == n, v.size(), A.rows());
    A.conservativeResize(n, m + 1);
    //A.topRightCorner(n, m) = A.topLeftCorner(n, m); // bad error, reason unknown!
    for (auto i = A.cols() - 1; i > 0; i--)
      A.col(i) = A.col(i - 1);
    A.col(0) = v;
  }
}

void removeColumnFromMatrix(
    Eigen::MatrixXd &A, int col)
{
  PRECICE_ASSERT(col < A.cols() && col >= 0, col, A.cols());
  for (int j = col; j < A.cols() - 1; j++)
    A.col(j) = A.col(j + 1);

  A.conservativeResize(A.rows(), A.cols() - 1);
}

void append(
    Eigen::VectorXd &v,
    double           value)
{
  int n = v.size();
  v.conservativeResize(n + 1);
  v(n) = value;
}

} // namespace precice::utils
