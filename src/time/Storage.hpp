#pragma once

#include <Eigen/Core>
#include "logging/Logger.hpp"

namespace precice::time {

class Storage {
public:
  /**
   * @brief Stores data samples in time and provides corresponding convenience functions.
   */
  Storage();

  /**
   * @brief Initialize storage by storing given value at time 0.0 and 1.0.
   *
   * @param values initial value
   */
  void initialize(Eigen::VectorXd values);

  /**
   * @brief Store a value in at a specific time.
   *
   * @param time the time associated with the value
   * @param value stored value
   */
  void setValueAtTime(double time, Eigen::VectorXd value);

  /**
   * @brief Get maximum normalized dt that is stored in this Storage.
   *
   * @return the maximum normalized dt from this Storage
   */
  double maxStoredNormalizedDt();

  /**
   * @brief Returns the value at time following "before" contained in this Storage.
   *
   * The stored normalized dt is larger or equal than "before". If "before" is a normalized dt stored in this Storage, this function returns the value at "before"
   *
   * @param before a double, where we want to find a normalized dt that comes directly after this one
   * @return Eigen::VectorXd a value in this Storage at or directly after "before"
   */
  Eigen::VectorXd getValueAtTimeAfter(double before);

  /**
   * @brief Get all normalized dts stored in this Storage sorted ascending.
   *
   * @return Eigen::VectorXd containing all stored normalized dts in ascending order.
   */
  Eigen::VectorXd getTimes();

  /**
   * @brief Get all normalized dts and values in ascending order (with respect to normalized dts)
   *
   * @return std::pair<Eigen::VectorXd, Eigen::MatrixXd> containing all stored times and values in ascending order (with respect to normalized dts).
   */
  std::pair<Eigen::VectorXd, Eigen::MatrixXd> getTimesAndValues();

  /**
   * @brief Number of stored times
   *
   * @return int number of stored times
   */
  int nTimes();

  /**
   * @brief Number of Dofs for each values
   *
   * @return int number of dofs
   */
  int nDofs();

  /**
   * @brief Move this Storage by storing the value at the end of the Storage at 0.0 and clearing the storage. Time 1.0 is initialized as value at 0.0
   */
  void move();

  /**
   * @brief Clear this Storage by deleting all values. If keepZero is true, keep values associated with 0.0.
   *
   * @param keepZero if true, keep value associated with 0.0.
   */
  void clear(bool keepZero = false);

private:
  /** @TODO Idea for more efficient data structure and redesign (do this when functionality is working and tested!)
   *   1. use Eigen::MatrixXd instead of map for _timeStepsStorage.
   *   2. create a member std::map<double, int> _timeSteps where (unique) time is mapped to column index of _timeStepsStorage that holds the corresponding sample. (Alternative: Use another Eigen::VectorXd to store times, but this enforces maintaining a consistent order for _timeSteps and _timeStepsStorage. This sounds complicated.)
   */
  /// Stores values on the current window associated with normalized dt.
  std::vector<std::pair<double, Eigen::VectorXd>> _sampleStorage;

  mutable logging::Logger _log{"time::Storage"};
};

} // namespace precice::time
