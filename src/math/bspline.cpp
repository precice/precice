#include "math/bspline.hpp"
#include "math/differences.hpp"
#include "utils/assertion.hpp"

#include <Eigen/Core>
#include <algorithm>
#include <cstdlib>
#include <functional>
#include <iterator>
#include <unsupported/Eigen/Splines>

namespace precice::math::bspline {

Eigen::VectorXd interpolateAt(Eigen::VectorXd ts, const Eigen::MatrixXd &xs, int splineDegree, double t)
{
  PRECICE_ASSERT(ts.size() >= 2, "Interpolation requires at least 2 samples");

  PRECICE_ASSERT(std::inner_product(
                     ts.data(), std::next(ts.data(), ts.size() - 1),
                     std::next(ts.data()), true,
                     std::logical_and<bool>{},
                     [](double l, double r) { return math::smaller(l, r); }),
                 "Timestamps must be strictly increasting", ts);

  // organize data in columns. Each column represents one sample in time.
  PRECICE_ASSERT(xs.cols() == ts.size());
  const int ndofs = xs.rows(); // number of dofs. Each dof needs it's own interpolant.

  // transform all times to the relative interval [0; 1]
  PRECICE_ASSERT(math::greaterEquals(t, ts(0)), "t is before the first sample!", t, ts(0), ts);                       // Only allowed to use BSpline for interpolation, not extrapolation.
  PRECICE_ASSERT(math::smallerEquals(t, ts(ts.size() - 1)), "t is after the last sample!", t, ts(ts.size() - 1), ts); // Only allowed to use BSpline for interpolation, not extrapolation.
  auto         relativeTime = [tsMin = ts(0), tsMax = ts(ts.size() - 1)](double t) -> double { return (t - tsMin) / (tsMax - tsMin); };
  const double tRelative    = relativeTime(t);
  ts                        = ts.unaryExpr(relativeTime);

  // interpolate all DOFs
  Eigen::VectorXd interpolated(ndofs);

  const int splineDimension = 1;

  // @todo implement cache to avoid unnecessary recomputation. Important! Need to reset cache when entering next window or iteration.
  for (int i = 0; i < ndofs; i++) {
    const auto spline = Eigen::SplineFitting<Eigen::Spline<double, splineDimension>>::Interpolate(xs.row(i), splineDegree, ts);
    interpolated[i]   = spline(tRelative)[0]; // get component of spline associated with xs.row(i)
  }

  return interpolated;
}

} // namespace precice::math::bspline
