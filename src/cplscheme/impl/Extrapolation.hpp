#pragma once

#include <Eigen/Core>
#include "logging/Logger.hpp"

namespace precice {

namespace testing {
// Forward declaration to friend the boost test struct
class ExtrapolationFixture;
} // namespace testing

namespace cplscheme {
namespace impl {

class Extrapolation {
  friend class testing::ExtrapolationFixture; // Make the fixture friend of this class
public:
  /**
   * @brief Extrapolation object which stores values of current and past time windows for performing extrapolation. 
   *
   * Storage still needs to be initialized with Extrapolation::initialize, before the Extrapolation can be used.
   *
   * @param extrapolationOrder defines the maximum extrapolation order supported by this extrapolation and reserves storage correspondingly.
   */
  Extrapolation(const int extrapolationOrder);

  /**
   * @brief Used to initialize _timeWindowsStorage according to required size.
   * @param valuesSize defines how many values one sample in time consists of
   */
  void initialize(const int valuesSize);

  /**
   * @brief Updates entry in _timeWindowsStorage corresponding to this window with given values
   * @param values new sample for this time window
   */
  void store(const Eigen::VectorXd &values);

  /**
   * @brief Called, when moving to the next time window. All entries in _timeWindowsStorage are shifted. The new entry is initialized as the value from the last window (= constant extrapolation)
   */
  void moveToNextWindow();

  /**
   * @brief getter for values at the current time window.
   */
  const Eigen::VectorXd getInitialGuess();

private:
  /// Set by initialize. Used for consistency checks.
  bool _storageIsInitialized = false;

  /// Stores values for several time windows.
  Eigen::MatrixXd _timeWindowsStorage;

  /// extrapolation order for this extrapolation
  int _extrapolationOrder; // @todo make const! Possible, if extrapolation order is set at configuration of data.

  /// number of stored samples in _timeWindowsStorage
  int _numberOfStoredSamples = 0;

  mutable logging::Logger _log{"cplscheme::Extrapolation"};

  /**
   * @brief returns number samples in time this extrapolation can store
   */
  int sizeOfSampleStorage();

  /**
   * @brief returns number of values per sample in time stored by this extrapolation
   */
  int valuesSize();

  /**
   * @brief Extrapolates values from _timeWindowsStorage using an extrapolation scheme of given order. 
   * 
   * If the order condition cannot be satisfied, since there are not enough samples available, the order is automatically reduced.
   * If order two is required, but only two samples are available, the extrapolation order is automatically reduced to one.
   */
  Eigen::VectorXd extrapolate();
};

} // namespace impl
} // namespace cplscheme
} // namespace precice
