#include "time/Waveform.hpp"
#include <algorithm>
#include <unsupported/Eigen/Splines>
#include "cplscheme/CouplingScheme.hpp"
#include "logging/LogMacros.hpp"
#include "math/differences.hpp"
#include "mesh/Data.hpp"
#include "time/Time.hpp"
#include "utils/EigenHelperFunctions.hpp"

namespace precice::time {

Waveform::Waveform(const int degree)
    : _degree(degree)
{
  PRECICE_ASSERT(Time::MIN_WAVEFORM_DEGREE <= _degree && _degree <= Time::MAX_WAVEFORM_DEGREE);
}

int Waveform::getDegree() const
{
  return _degree;
}

time::Storage &Waveform::timeStepsStorage()
{
  return _timeStepsStorage;
}

const time::Storage &Waveform::timeStepsStorage() const
{
  return _timeStepsStorage;
}

// helper function to compute x(t) from given data (x0,t0), (x1,t1), ..., (xn,tn) via B-spline interpolation (implemented using Eigen).
Eigen::VectorXd bSplineInterpolationAt(double t, Eigen::VectorXd ts, Eigen::MatrixXd xs, int splineDegree)
{
  // organize data in columns. Each column represents one sample in time.
  PRECICE_ASSERT(xs.cols() == ts.size());
  const int ndofs = xs.rows(); // number of dofs. Each dof needs it's own interpolant.

  double tRelative = (t - ts(0)) / (ts(ts.size() - 1) - ts(0)); // Eigen requires us to scale the time t to the interval [0,1]

  Eigen::VectorXd tsRelative(ts.size());

  for (int i = 0; i < ts.size(); i++) {
    tsRelative[i] = (ts(i) - ts(0)) / (ts(ts.size() - 1) - ts(0)); // Eigen requires us to scale the knot vector ts to the interval [0,1]
  }

  Eigen::VectorXd interpolated(ndofs);

  const int splineDimension = 1;

  // @todo implement cache to avoid unnecessary recomputation. Important! Need to reset cache when entering next window or iteration.
  for (int i = 0; i < ndofs; i++) {
    const auto spline = Eigen::SplineFitting<Eigen::Spline<double, splineDimension>>::Interpolate(xs.row(i), splineDegree, tsRelative);
    interpolated[i]   = spline(tRelative)[0]; // get component of spline associated with xs.row(i)
  }

  return interpolated;
}

Eigen::VectorXd Waveform::sample(double time) const
{
  const int usedDegree = computeUsedDegree(_degree, _timeStepsStorage.nTimes());

  if (usedDegree == 0) {
    return this->_timeStepsStorage.getValuesAtOrAfter(time);
  }

  const auto data = _timeStepsStorage.getTimesAndValues();

  return bSplineInterpolationAt(time, data.first, data.second, usedDegree);
}

int Waveform::computeUsedDegree(int requestedDegree, int numberOfAvailableSamples) const
{
  int usedDegree = -1;
  PRECICE_ASSERT(requestedDegree <= 3);
  if (requestedDegree == 0 || numberOfAvailableSamples < 2) {
    usedDegree = 0;
  } else if (requestedDegree == 1 || numberOfAvailableSamples < 3) {
    usedDegree = 1;
  } else if (requestedDegree == 2 || numberOfAvailableSamples < 4) {
    usedDegree = 2;
  } else if (requestedDegree == 3 || numberOfAvailableSamples < 5) {
    usedDegree = 3;
  } else {
    PRECICE_ASSERT(false); // not supported
  }
  return usedDegree;
}

} // namespace precice::time
