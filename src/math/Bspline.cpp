#include "math/Bspline.hpp"
#include "math/differences.hpp"
#include "utils/assertion.hpp"

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <algorithm>
#include <cstdlib>
#include <unsupported/Eigen/Splines>
#include "math/differences.hpp"
#include "utils/assertion.hpp"

namespace precice::math {

Bspline::Bspline(Eigen::VectorXd ts, const Eigen::MatrixXd &xs, int splineDegree)
{

  PRECICE_ASSERT(ts.size() >= 2, "Interpolation requires at least 2 samples");
  PRECICE_ASSERT(std::is_sorted(ts.begin(), ts.end()), "Timestamps must be sorted");

  // organize data in columns. Each column represents one sample in time.
  PRECICE_ASSERT(xs.cols() == ts.size());
  _ndofs            = xs.rows(); // number of dofs. Each dof needs its own interpolant.
  _tsMin            = ts(0);
  _tsMax            = ts(ts.size() - 1);
  auto relativeTime = [tsMin = _tsMin, tsMax = _tsMax](double t) -> double { return (t - tsMin) / (tsMax - tsMin); };
  ts                = ts.unaryExpr(relativeTime);

  //The code for computing the knots and the control points is copied from Eigens bspline interpolation with minor modifications, https://gitlab.com/libeigen/eigen/-/blob/master/unsupported/Eigen/src/Splines/SplineFitting.h

  Eigen::KnotAveraging(ts, splineDegree, _knots);
  Eigen::DenseIndex                   n = xs.cols();
  std::vector<Eigen::Triplet<double>> matrixEntries;
  matrixEntries.reserve(n * splineDegree + 2);

  for (Eigen::DenseIndex i = 1; i < n - 1; ++i) {
    const Eigen::DenseIndex span      = Eigen::Spline<double, 1>::Span(ts[i], splineDegree, _knots);
    auto                    basisFunc = Eigen::Spline<double, 1>::BasisFunctions(ts[i], splineDegree, _knots);

    for (Eigen::DenseIndex j = 0; j < splineDegree + 1; ++j) {
      matrixEntries.emplace_back(i, span - splineDegree + j, basisFunc(j));
    }
  }

  matrixEntries.emplace_back(0, 0, 1.0);
  matrixEntries.emplace_back(n - 1, n - 1, 1.0);

  Eigen::SparseMatrix<double> A(n, n);
  A.setFromTriplets(matrixEntries.begin(), matrixEntries.end());
  A.makeCompressed();

  Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> qr;
  qr.analyzePattern(A);
  qr.factorize(A);

  _ctrls = qr.solve(xs.transpose());
}

Eigen::VectorXd Bspline::interpolateAt(double t) const
{
  // transform t to the relative interval [0; 1]
  const double tRelative = std::clamp((t - _tsMin) / (_tsMax - _tsMin), 0.0, 1.0);

  Eigen::VectorXd interpolated(_ndofs);
  constexpr int   splineDimension = 1;

  for (int i = 0; i < _ndofs; i++) {
    interpolated[i] = Eigen::Spline<double, splineDimension>(_knots, _ctrls.col(i))(tRelative)[0];
  }

  return interpolated;
}
} // namespace precice::math
