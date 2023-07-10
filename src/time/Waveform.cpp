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

Waveform::Waveform(
    const int degree, mesh::PtrData data)
    : _degree(degree), _data(data)
{
  PRECICE_ASSERT(Time::MIN_WAVEFORM_DEGREE <= _degree && _degree <= Time::MAX_WAVEFORM_DEGREE);
}

int Waveform::getDegree() const
{
  return _degree;
}

void Waveform::setDegree(int degree)
{
  _degree = degree;
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
    const auto spline = Eigen::SplineFitting<Eigen::Spline<double, splineDimension>>::Interpolate(xs.row(i), splineDegree, ts);
    interpolated[i]   = spline(t)[0]; // get component of spline associated with xs.row(i)
  }

  return interpolated;
}

Eigen::VectorXd Waveform::sample(double normalizedDt) const
{
  const int usedDegree = computeUsedDegree(_degree, _data->timeStepsStorage().nTimes());

  PRECICE_ASSERT(math::equals(this->_data->timeStepsStorage().maxStoredNormalizedDt(), time::Storage::WINDOW_END), this->_data->timeStepsStorage().maxStoredNormalizedDt()); // sampling is only allowed, if a window is complete.

  if (_degree == 0) {
    return this->_data->timeStepsStorage().getValuesAtOrAfter(normalizedDt);
  }

  PRECICE_ASSERT(usedDegree >= 1);

  const auto data = _data->timeStepsStorage().getTimesAndValues();

  return bSplineInterpolationAt(normalizedDt, data.first, data.second, usedDegree);
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
