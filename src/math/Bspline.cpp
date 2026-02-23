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
#include <chrono>  // Add for timing
#include <iostream> // Add for debug output (optional)

namespace precice::math {

// Add this struct at the top or inside the cpp file for collecting metrics
struct BsplinePerformanceMetrics {
    double constructionTime = 0.0;
    double factorizationTime = 0.0;
    double solveTime = 0.0;
    double totalTime = 0.0;
    size_t matrixSize = 0;
    int splineDegree = 0;
    
    void print() const {
        std::cout << "Matrix " << matrixSize << "x" << matrixSize 
                  << " (degree " << splineDegree << "): "
                  << "construct=" << constructionTime << "ms, "
                  << "factor=" << factorizationTime << "ms, "
                  << "solve=" << solveTime << "ms, "
                  << "total=" << totalTime << "ms" << std::endl;
    }
};

// Optional: Enable/disable profiling via environment variable
static bool enableProfiling() {
    static bool enabled = []() {
        const char* env = std::getenv("PRECICE_BSPLINE_PROFILE");
        return env && std::string(env) == "1";
    }();
    return enabled;
}

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

  // The code for computing the knots and the control points is copied from Eigens bspline interpolation with some modifications
  // https://gitlab.com/libeigen/eigen/-/blob/master/unsupported/Eigen/src/Splines/SplineFitting.h

  // Start timing if profiling enabled
  auto totalStart = std::chrono::high_resolution_clock::now();
  auto knotStart = totalStart;

  // 1. Compute the knot vector
  Eigen::KnotAveraging(ts, splineDegree, _knots);
  
  auto knotEnd = std::chrono::high_resolution_clock::now();
  double knotTime = std::chrono::duration<double, std::milli>(knotEnd - knotStart).count();

  // 2. Compute the control points
  // We use a nxn sparse matrix with 2 + (n-2) * (d+1) entries and thus a fill-factor < 0.5.
  Eigen::DenseIndex                   n = xs.cols();
  std::vector<Eigen::Triplet<double>> matrixEntries;
  matrixEntries.reserve(2 + (n - 2) * (splineDegree + 1));

  auto constructStart = std::chrono::high_resolution_clock::now();

  matrixEntries.emplace_back(0, 0, 1.0);
  for (Eigen::DenseIndex i = 1; i < n - 1; ++i) {
    const Eigen::DenseIndex span      = Eigen::Spline<double, 1>::Span(ts[i], splineDegree, _knots);
    auto                    basisFunc = Eigen::Spline<double, 1>::BasisFunctions(ts[i], splineDegree, _knots);

    for (Eigen::DenseIndex j = 0; j < splineDegree + 1; ++j) {
      matrixEntries.emplace_back(i, span - splineDegree + j, basisFunc(j));
    }
  }
  matrixEntries.emplace_back(n - 1, n - 1, 1.0);
  PRECICE_ASSERT(matrixEntries.capacity() == matrixEntries.size(), matrixEntries.capacity(), matrixEntries.size(), n, splineDegree);

  Eigen::SparseMatrix<double> A(n, n);
  A.setFromTriplets(matrixEntries.begin(), matrixEntries.end());
  A.makeCompressed();

  auto constructEnd = std::chrono::high_resolution_clock::now();
  double constructTime = std::chrono::duration<double, std::milli>(constructEnd - constructStart).count();

  Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> qr;
  
  auto factorStart = std::chrono::high_resolution_clock::now();
  qr.analyzePattern(A);
  qr.factorize(A);
  auto factorEnd = std::chrono::high_resolution_clock::now();
  double factorTime = std::chrono::duration<double, std::milli>(factorEnd - factorStart).count();

  auto solveStart = std::chrono::high_resolution_clock::now();
  _ctrls = qr.solve(xs.transpose());
  auto solveEnd = std::chrono::high_resolution_clock::now();
  double solveTime = std::chrono::duration<double, std::milli>(solveEnd - solveStart).count();

  auto totalEnd = std::chrono::high_resolution_clock::now();
  double totalTime = std::chrono::duration<double, std::milli>(totalEnd - totalStart).count();

  // Print metrics if profiling enabled
  if (enableProfiling()) {
    BsplinePerformanceMetrics metrics;
    metrics.constructionTime = constructTime;
    metrics.factorizationTime = factorTime;
    metrics.solveTime = solveTime;
    metrics.totalTime = totalTime;
    metrics.matrixSize = n;
    metrics.splineDegree = splineDegree;
    metrics.print();
  }

  // Optional: Add fill factor analysis for debugging
  if (enableProfiling() && n < 100) { // Only for small matrices to avoid spam
    double nonZeros = 2 + (n - 2) * (splineDegree + 1);
    double fillFactor = nonZeros / (n * n);
    std::cout << "  Fill factor: " << fillFactor * 100 << "% (" 
              << nonZeros << "/" << n*n << " non-zeros)" << std::endl;
  }
}

Eigen::VectorXd Bspline::interpolateAt(double t) const
{
  // transform t to the relative interval [0; 1]
  const double tRelative = std::clamp((t - _tsMin) / (_tsMax - _tsMin), 0.0, 1.0);

  Eigen::VectorXd interpolated(_ndofs);
  constexpr int   splineDimension = 1;

  // Optional: Time interpolation if profiling enabled
  auto start = std::chrono::high_resolution_clock::now();

  for (int i = 0; i < _ndofs; i++) {
    interpolated[i] = Eigen::Spline<double, splineDimension>(_knots, _ctrls.col(i))(tRelative)[0];
  }

  if (enableProfiling()) {
    auto end = std::chrono::high_resolution_clock::now();
    double time = std::chrono::duration<double, std::milli>(end - start).count();
    std::cout << "Interpolation at t=" << t << " took " << time << "ms" << std::endl;
  }

  return interpolated;
}

} // namespace precice::math