#pragma once

#include <Eigen/Core>
#include "logging/Logger.hpp"
#include "time/Stample.hpp"

namespace precice::time {

class Storage {
public:
  /// Fixed time associated with beginning of window
  static const double WINDOW_START;

  /// Fixed time associated with end of window
  static const double WINDOW_END;

  /**
   * @brief Stores data samples in time and provides corresponding convenience functions.
   *
   * The Storage must be initialized before it can be used. Then values can be stored in the Storage. It is only allowed to store samples with increasing times. Overwriting existing samples or writing samples with a time smaller then the maximum stored time is forbidden.
   * The Storage is considered complete, when a sample with time 1.0 is provided. Then one can only sample from the storage. To add further samples one needs to trim or clear the storage first.
   *
   * This Storage is used in the context of Waveform relaxation where samples in time are provided. Starting at the beginning of the window with time 0.0 and reaching the end of the window with time 1.0.
   */
  Storage();

  /**
   * @brief Initialize storage by storing given sample at time 0.0 and 1.0.
   *
   * @param sample initial sample
   */
  void initialize(time::Sample sample);

  /**
   * @brief Store Sample at a specific time.
   *
   * It is only allowed to store a Sample in time that comes after a Sample that was already stored. Therefore, time has to be larger than maxStoredNormalizedDt. Overwriting existing samples is forbidden. The function trim() or clear() should be used to be able provide new samples.
   *
   * @param time the time associated with the sample
   * @param sample stored sample
   */
  void setSampleAtTime(double time, Sample sample);

  /**
   * @brief Get maximum normalized dt that is stored in this Storage.
   *
   * @return the maximum normalized dt from this Storage
   */
  double maxStoredNormalizedDt() const;

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
   * @brief Get the vector of Stamples
   *
   * @return std::vector<Stample>
   */
  const std::vector<Stample> &getStamples() const;

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
   * @brief Move this Storage by storing the values at the end of the Storage at 0.0 and clearing the storage. Time 1.0 is initialized as values at 0.0
   */
  void move();

  /**
   * @brief Trims this Storage by deleting all values except values associated with 0.0.
   */
  void trim();

  /**
   * @brief Clear this Storage by deleting all values.
   */
  void clear();

private:
  /// Stores Stamples on the current window
  std::vector<Stample> _stampleStorage;

  mutable logging::Logger _log{"time::Storage"};
};

} // namespace precice::time
