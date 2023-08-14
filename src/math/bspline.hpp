#pragma once

#include <Eigen/Core>

namespace precice::math::bspline {

/** Uses B-Spline interpolation to compute compute x(t) from given data (x0,t0), (x1,t1), ..., (xn,tn).
 *
 * @pre Requires at least 2 samples.
 * @pre Timestamps must be sorted from lowest to highest.
 * @pre splineDegree must be larger than 0.
 * @pre t must be within [ts.min; ts.max].
 */
Eigen::VectorXd interpolateAt(Eigen::VectorXd ts, const Eigen::MatrixXd &xs, int splineDegree, double t);

} // namespace precice::math::bspline
