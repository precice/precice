#pragma once
#include <Eigen/Core>

namespace precice::math {

class Bspline {

public:
  /**
 * @brief Initialises the B-Spline interpolation with the given data (x0,t0), (x1,t1), ..., (xn,tn) and computes the knots and control points.
 * The code for computing the knots and the control points is copied from Eigens bspline interpolation with minor modifications, https://gitlab.com/libeigen/eigen/-/blob/master/unsupported/Eigen/src/Splines/SplineFitting.h
 *
 * @param ts the timestamps which must be sorted from lowest to highest and contain at least 2 samples.
 * @param xs the data to be interpolated. It has to contain at least two samples.
 * @param splineDegree the used spline degree, which has to be larger than 0
 */
  Bspline(Eigen::VectorXd ts, const Eigen::MatrixXd &xs, int splineDegree);

  /**
 * @brief Samples the B-Spline interpolation
 *
 * @param t must be within [_tsMin; _tsMax].
 * @return the inteprolant x(t) if t does not equal any of the timestamps in ts.
 * If t equals any of the timestamps in ts the corresponding value to t is returned.
 */

  Eigen::VectorXd interpolateAt(double t) const;

private:
  Eigen::VectorXd _knots; // Cache to store previously computed knots
  Eigen::MatrixXd _ctrls; // Cache to store previously computed control points
  double          _tsMin; // The minimal time of the bspline
  double          _tsMax; // The maximal time of the bspline
  int             _ndofs; // The degrees of freedom of the data
};
} // namespace precice::math
