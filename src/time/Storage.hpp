#pragma once

#include <Eigen/Core>
#include "logging/Logger.hpp"

namespace precice {
namespace time {

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
   * @brief Get the Value stored at time
   *
   * @param time time for which we are requesting a stored value
   * @return Eigen::VectorXd value associated with time
   */
  Eigen::VectorXd getValueAtTime(double time);

  /**
   * @brief Store a value in at a specific time.
   *
   * @param time the time associated with the value
   * @param value stored value
   */
  void setValueAtTime(double time, Eigen::VectorXd value);

  /**
   * @brief Get maximum time is stored in this Storage.
   *
   * @return the maximum time from found in this Storage
   */
  double maxTime();

  /**
   * @brief Returns the stored time closest to "before" contained in this Storage.
   *
   * The stored time is larger or equal than "before". If "before" is a time stored in this Storage, this function returns "before"
   *
   * @param before a double, where we want to find a time that comes directly after this one
   * @return double a time in this Storage which is larger or equal to "before"
   */
  double getClosestTimeAfter(double before);

  /**
   * @brief Get all times stored in this Storage sorted ascending.
   *
   * @return Eigen::VectorXd containing all stored times in ascending order.
   */
  Eigen::VectorXd getTimes();

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
   * @brief Clear this Storage by deleting all values. If keepZero is true, keeps Value associated with 0.0.
   *
   * @param keepZero if true, keep value associated with 0.0.
   */
  void clear(bool keepZero = false);

private:
  /** @TODO Idea for more efficient data structure and redesign (do this when functionality is working and tested!)
   *   1. use Eigen::MatrixXd instead of map for _timeStepsStorage.
   *   2. create a member std::map<double, int> _timeSteps where (unique) time is mapped to column index of _timeStepsStorage that holds the corresponding sample. (Alternative: Use another Eigen::VectorXd to store times, but this enforces maintaining a consistent order for _timeSteps and _timeStepsStorage. This sounds complicated.)
   */
  /// Stores values on the current window.
  std::map<double, Eigen::VectorXd> _storageDict;

  mutable logging::Logger _log{"time::Storage"};
};

} // namespace time
} // namespace precice
