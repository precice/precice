#pragma once

#include <Eigen/Core>
#include <boost/range.hpp>
#include <optional>
#include "logging/Logger.hpp"
#include "math/Bspline.hpp"
#include "time/Stample.hpp"

namespace precice::time {

class Storage {
public:
  /**
   * @brief Stores data samples in time and provides corresponding convenience functions.
   *
   * The Storage must be initialized before it can be used. Then values can be stored in the Storage. It is only allowed to store samples with increasing times. Overwriting existing samples or writing samples with a time smaller then the maximum stored time is forbidden.
   * The Storage is considered complete, when a sample for the end of the current window is provided. Then one can only sample from the storage. To add further samples one needs to trim the storage or move to the next time window first.
   *
   * This Storage is used in the context of Waveform relaxation where samples in time are provided.
   */
  Storage();

  /**
   * @brief Copy assignment operator to assign Storage to this Storage
   *
   * @param other Storage
   * @return Storage&
   */
  Storage &operator=(const Storage &other);

  /**
   * @brief Store Sample at a specific time.
   *
   * It is only allowed to store a Sample in time that comes after a Sample that was already stored. Therefore, time has to be larger than maxStoredTime. Overwriting existing samples is forbidden. The function trim() should be used before providing new samples.
   *
   * @param time the time associated with the sample
   * @param sample stored sample
   */
  void setSampleAtTime(double time, Sample sample);

  void setInterpolationDegree(int interpolationDegree);

  int getInterpolationDegree() const;

  /**
   * @brief Get maximum time that is stored in this Storage.
   *
   * @return the maximum time from this Storage
   */
  double maxStoredTime() const;

  /**
   * @brief Returns the Sample at time following "before" contained in this Storage.
   *
   * The stored normalized dt is larger or equal than "before". If "before" is a normalized dt stored in this Storage, this function returns the Sample at "before"
   *
   * @param before a double, where we want to find a normalized dt that comes directly after this one
   * @return Sample in this Storage at or directly after "before"
   */
  Sample getSampleAtOrAfter(double before) const;

  /**
   * @brief Get all normalized dts stored in this Storage sorted ascending.
   *
   * @return Eigen::VectorXd containing all stored normalized dts in ascending order.
   */
  Eigen::VectorXd getTimes() const;

  /**
   * @brief Get the stamples
   *
   * @return boost range of stamples
   */
  auto stamples() const
  {
    return boost::make_iterator_range(_stampleStorage);
  }

  auto stamples()
  {
    return boost::make_iterator_range(_stampleStorage);
  }

  bool empty() const;

  const time::Stample &last() const;

  /**
   * @brief Get all normalized dts and values in ascending order (with respect to normalized dts)
   *
   * @return std::pair<Eigen::VectorXd, Eigen::MatrixXd> containing all stored times and values in ascending order (with respect to normalized dts).
   */
  std::pair<Eigen::VectorXd, Eigen::MatrixXd> getTimesAndValues() const;

  /**
   * @brief Number of stored times
   *
   * @return int number of stored times
   */
  int nTimes() const;

  /**
   * @brief Number of Dofs for each values
   *
   * @return int number of dofs
   */
  int nDofs() const;

  /**
   * @brief Move this Storage by deleting all stamples except the one at the end of the window.
   */
  void move();

  /**
   * @brief Trims this Storage by deleting all values except values associated with the window start.
   */
  void trim();

  /**
   * @brief Clears this Storage by deleting all values.
   */
  void clear();

  void trimBefore(double time);

  void trimAfter(double time);

  /**
   * @brief Need to use interpolation for the case with changing time grids
   *
   * @param time a double, where we want to sample the waveform
   * @return Eigen::VectorXd values in this Storage at or directly after "before"
  */
  Eigen::VectorXd sample(double time) const;

  Eigen::MatrixXd sampleGradients(double time) const;

private:
  /// Stores Stamples on the current window
  std::vector<Stample> _stampleStorage;

  mutable logging::Logger _log{"time::Storage"};

  int _degree;

  mutable std::optional<math::Bspline> _bspline;

  /**
   * @brief Computes which degree may be used for interpolation.
   *
   * Actual degree of interpolating B-spline is determined by number of stored samples and maximum degree defined by the user.
   * Example: If only two samples are available, the maximum degree we may use is 1, even if the user demands degree 2.
   *
   * @param requestedDegree B-spline degree requested by the user.
   * @param numberOfAvailableSamples Samples available for interpolation.
   * @return B-spline degree that may be used.
   */
  int computeUsedDegree(int requestedDegree, int numberOfAvailableSamples) const;

  time::Sample getSampleAtBeginning();

  time::Sample getSampleAtEnd();

  int findTimeId(double time) const;
};

} // namespace precice::time
