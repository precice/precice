#pragma once

#include <Eigen/Core>
#include <boost/range.hpp>
#include "logging/Logger.hpp"
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
   * @brief Initialize storage by storing given sample at window start and window end.
   *
   * @param sample initial sample
   */
  void initialize(time::Sample sample);

  /**
   * @brief Store Sample at a specific time.
   *
   * It is only allowed to store a Sample in time that comes after a Sample that was already stored. Therefore, time has to be larger than maxStoredTime. Overwriting existing samples is forbidden. The function trim() should be used before providing new samples.
   *
   * @param time the time associated with the sample
   * @param sample stored sample
   */
  void setSampleAtTime(double time, Sample sample);

  /**
   * @brief Get maximum time that is stored in this Storage.
   *
   * @return the maximum time from this Storage
   */
  double maxStoredTime() const;

  /**
   * @brief Returns the values at time following "before" contained in this Storage.
   *
   * The stored normalized dt is larger or equal than "before". If "before" is a normalized dt stored in this Storage, this function returns the values at "before"
   *
   * @param before a double, where we want to find a normalized dt that comes directly after this one
   * @return Eigen::VectorXd values in this Storage at or directly after "before"
   */
  Eigen::VectorXd getValuesAtOrAfter(double before) const;

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

private:
  /// Stores Stamples on the current window
  std::vector<Stample> _stampleStorage;

  mutable logging::Logger _log{"time::Storage"};

  /// End time of the current window
  double _currentWindowStart;

  time::Sample getSampleAtBeginning();

  time::Sample getSampleAtEnd();
};

} // namespace precice::time
