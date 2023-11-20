#include "math/Bspline.hpp"
#include "math/differences.hpp"
#include "utils/assertion.hpp"

#include <Eigen/Core>
#include <algorithm>
#include <cstdlib>
#include <unsupported/Eigen/Splines>
#include "math/differences.hpp"
#include "utils/assertion.hpp"

namespace precice::math {

Bspline::Bspline(Eigen::VectorXd ts, const Eigen::MatrixXd &xs, int splineDegree)
{

  PRECICE_ASSERT(ts.size() >= 2, "Interpolation requires at least 2 samples");
#if EIGEN_VERSION_AT_LEAST(3, 4, 0)
  PRECICE_ASSERT(std::is_sorted(ts.begin(), ts.end()), "Timestamps must be sorted");
#else
  PRECICE_ASSERT(std::is_sorted(ts.data(), ts.data() + ts.size()), "Timestamps must be sorted");
#endif

  // organize data in columns. Each column represents one sample in time.
  PRECICE_ASSERT(xs.cols() == ts.size());
  _ndofs            = xs.rows(); // number of dofs. Each dof needs its own interpolant.
  _tsMin            = ts(0);
  _tsMax            = ts(ts.size() - 1);
  auto relativeTime = [tsMin = _tsMin, tsMax = _tsMax](double t) -> double { return (t - tsMin) / (tsMax - tsMin); };
  ts                = ts.unaryExpr(relativeTime);

  //The code for computing the knots and the control points is copied from Eigens bspline interpolation with minor modifications, https://gitlab.com/libeigen/eigen/-/blob/master/unsupported/Eigen/src/Splines/SplineFitting.h

  Eigen::KnotAveraging(ts, splineDegree, _knots);
  Eigen::DenseIndex n = xs.cols();
  Eigen::MatrixXd   A = Eigen::MatrixXd::Zero(n, n);

  for (Eigen::DenseIndex i = 1; i < n - 1; ++i) {
    //Attempt at hack... the spline dimension is not used here explicitly
    const Eigen::DenseIndex span = Eigen::Spline<double, 1>::Span(ts[i], splineDegree, _knots);

    // The segment call should somehow be told the spline order at compile time.
    A.row(i).segment(span - splineDegree, splineDegree + 1) = Eigen::Spline<double, 1>::BasisFunctions(ts[i], splineDegree, _knots);
  }
  A(0, 0)         = 1.0;
  A(n - 1, n - 1) = 1.0;

  auto qr = A.householderQr();

  Eigen::MatrixXd controls = Eigen::MatrixXd::Zero(n, _ndofs);

  for (int i = 0; i < _ndofs; i++) {
    controls.col(i) = qr.solve(xs.row(i).transpose());
  }
  _ctrls = std::move(controls);
}

Eigen::VectorXd Bspline::interpolateAt(double t) const
{
  if (math::equals(t, _tsMax)) {
    t = _tsMax; // set t exactly to _tsMax, to avoid artifacts originating from floating-point arithmetic
  } else if (math::equals(t, _tsMin)) {
    t = _tsMin; // set t exactly to _tsMin, to avoid artifacts originating from floating-point arithmetic
  }

  // transform t to the relative interval [0; 1]
  const double tRelative = (t - _tsMin) / (_tsMax - _tsMin);
  PRECICE_ASSERT(math::greaterEquals(tRelative, 0.0), "t is before the first sample!", tRelative, 0.0, _tsMin); // Only allowed to use BSpline for interpolation, not extrapolation.
  PRECICE_ASSERT(math::smallerEquals(tRelative, 1.0), "t is after the last sample!", tRelative, 1.0, _tsMax);   // Only allowed to use BSpline for interpolation, not extrapolation.

  Eigen::VectorXd interpolated(_ndofs);
  constexpr int   splineDimension = 1;

  for (int i = 0; i < _ndofs; i++) {
    interpolated[i] = Eigen::Spline<double, splineDimension>(_knots, _ctrls.col(i))(tRelative)[0];
  }

  return interpolated;
}
} // namespace precice::math
