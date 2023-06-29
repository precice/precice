#include "time/Waveform.hpp"
#include <algorithm>
#include "cplscheme/CouplingScheme.hpp"
#include "logging/LogMacros.hpp"
#include "math/differences.hpp"
#include "mesh/Data.hpp"
#include "time/Time.hpp"
#include "utils/EigenHelperFunctions.hpp"
#include "utils/assertion.hpp"

#include <cstddef>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>

namespace {
// Fit B-spline to positions and values and interpolate at a given position
double fitAndInterpolateBSpline(const Eigen::VectorXd &positions, const Eigen::VectorXd &values, int degree, double sampleAt)
{
  // Check if the input vectors have the same size
  PRECICE_ASSERT(positions.size() > 1, "There are 2 samples required");
  PRECICE_ASSERT(positions.size() == values.size(), "Input vectors must have the same size.");
  PRECICE_ASSERT(degree >= 0);
  PRECICE_ASSERT(degree <= positions.size(), "Cannot fit if degree > nsamples")

  const size_t minBreakpoints = positions.size() + 2 - degree;

  // Initialize the B-spline workspace
  const size_t           numPositions   = positions.size();
  const size_t           numBreakpoints = std::min<size_t>(13, minBreakpoints); // will always be 2
  gsl_bspline_workspace *workspace      = gsl_bspline_alloc(degree, numBreakpoints);
  const size_t           numKnots       = gsl_bspline_ncoeffs(workspace);

  // Compute the knots
  gsl_bspline_knots_uniform(positions[0], positions[numPositions - 1], workspace);

  // Create the matrix for the B-spline fit
  gsl_matrix *X = gsl_matrix_alloc(numPositions, numKnots);
  for (size_t i = 0; i < numPositions; ++i) {
    double          x   = positions[i];
    gsl_vector_view row = gsl_matrix_row(X, i);
    gsl_bspline_eval(x, &row.vector, workspace);
  }

  // Create the vector for the values
  gsl_vector *Y = gsl_vector_alloc(numPositions);
  for (size_t i = 0; i < numPositions; ++i) {
    gsl_vector_set(Y, i, values[i]);
  }

  // Perform the B-spline fit
  gsl_vector *                   c            = gsl_vector_alloc(numKnots);
  gsl_matrix *                   cov          = gsl_matrix_alloc(numKnots, numKnots);
  gsl_multifit_linear_workspace *fitWorkspace = gsl_multifit_linear_alloc(numPositions, numKnots);
  [[maybe_unused]] double        chisq;
  gsl_multifit_linear(X, Y, c, cov, &chisq, fitWorkspace);

  // Interpolate the B-spline at the given position
  double      interpolatedValue = 0.0;
  gsl_vector *basis             = gsl_vector_alloc(numKnots);
  gsl_bspline_eval(sampleAt, basis, workspace);
  for (size_t i = 0; i < numKnots; ++i) {
    double basisValue  = gsl_vector_get(basis, i);
    double coefficient = gsl_vector_get(c, i);
    interpolatedValue += basisValue * coefficient;
  }

  // Clean up
  gsl_bspline_free(workspace);
  gsl_matrix_free(X);
  gsl_vector_free(Y);
  gsl_vector_free(c);
  gsl_matrix_free(cov);
  gsl_multifit_linear_free(fitWorkspace);
  gsl_vector_free(basis);

  return interpolatedValue;
}
} // namespace

namespace precice::time {

Waveform::Waveform(
    const int interpolationOrder, mesh::PtrData data)
    : _interpolationOrder(interpolationOrder), _data(data)
{
  PRECICE_ASSERT(Time::MIN_INTERPOLATION_ORDER <= _interpolationOrder && _interpolationOrder <= Time::MAX_INTERPOLATION_ORDER);
}

int Waveform::getInterpolationOrder() const
{
  return _interpolationOrder;
}

// helper function to compute x(t) from given data (x0,t0), (x1,t1), ..., (xn,tn) via B-spline interpolation (implemented using Eigen).
Eigen::VectorXd bSplineInterpolationAt(double t, Eigen::VectorXd ts, Eigen::MatrixXd xs, int splineDegree)
{
  // organize data in columns. Each column represents one sample in time.
  PRECICE_ASSERT(xs.cols() == ts.size());
  const int ndofs = xs.rows(); // number of dofs. Each dof needs it's own interpolant.

  Eigen::VectorXd interpolated(ndofs);

  const int splineDimension = 1;

  // @todo implement cache to avoid unnecessary recomputation. Important! Need to reset cache when entering next window or iteration.
  for (int i = 0; i < ndofs; i++) {
    interpolated[i] = fitAndInterpolateBSpline(ts, xs.row(i), splineDegree, t);
  }

  return interpolated;
}

Eigen::VectorXd Waveform::sample(double normalizedDt) const
{
  const int usedOrder = computeUsedOrder(_interpolationOrder, _data->timeStepsStorage().nTimes());

  PRECICE_ASSERT(math::equals(this->_data->timeStepsStorage().maxStoredNormalizedDt(), time::Storage::WINDOW_END), this->_data->timeStepsStorage().maxStoredNormalizedDt()); // sampling is only allowed, if a window is complete.

  if (_interpolationOrder == 0) {
    return this->_data->timeStepsStorage().getValuesAtOrAfter(normalizedDt);
  }

  PRECICE_ASSERT(usedOrder >= 1);

  const auto data = _data->timeStepsStorage().getTimesAndValues();

  return bSplineInterpolationAt(normalizedDt, data.first, data.second, usedOrder);
}

int Waveform::computeUsedOrder(int requestedOrder, int numberOfAvailableSamples) const
{
  int usedOrder = -1;
  PRECICE_ASSERT(requestedOrder <= 3);
  if (requestedOrder == 0 || numberOfAvailableSamples < 2) {
    usedOrder = 0;
  } else if (requestedOrder == 1 || numberOfAvailableSamples < 3) {
    usedOrder = 1;
  } else if (requestedOrder == 2 || numberOfAvailableSamples < 4) {
    usedOrder = 2;
  } else if (requestedOrder == 3 || numberOfAvailableSamples < 5) {
    usedOrder = 3;
  } else {
    PRECICE_ASSERT(false); // not supported
  }
  return usedOrder;
}

} // namespace precice::time
