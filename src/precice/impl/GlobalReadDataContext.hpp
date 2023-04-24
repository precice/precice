#pragma once

#include <Eigen/Core>

#include "DataContext.hpp"
#include "logging/Logger.hpp"
#include "time/SharedPointer.hpp"
#include "time/Time.hpp"

namespace precice {
namespace impl {

/**
 * @brief Stores one Global (meshless) Data object. Context stores data to be read from. Additionally stores Waveform object associated with _providedData.
 *
 * Derived from DataContext
 */
class GlobalReadDataContext : public DataContext {
public:
  /**
   * @brief Construct a new GlobalReadDataContext object without a mapping.
   *
   * @param data Data associated with this GlobalReadDataContext.
   * @param interpolationOrder Order of the Waveform stored by this GlobalReadDataContext.
   */
  GlobalReadDataContext(
      mesh::PtrGlobalData data,
      int                 interpolationOrder = time::Time::DEFAULT_INTERPOLATION_ORDER);

  /**
   * @brief Gets _interpolationOrder of _waveform
   *
   * @return _interpolationOrder of _waveform
   */
  int getInterpolationOrder() const;

  /**
   * @brief Samples data at a given point in time within the current time window
   *
   * @param normalizedDt Point in time where waveform is sampled. Must be normalized to [0,1], where 0 refers to the beginning and 1 to the end of the current time window.
   */
  Eigen::VectorXd sampleWaveformAt(double normalizedDt);

  /**
   * @brief Initializes the _waveform as a constant function with values from _providedData.
   */
  void initializeWaveform();

  /**
   * @brief Updates _waveform when moving to the next time window.
   */
  void moveToNextWindow();

  /**
   * @brief Stores _providedData as first sample of _waveform.
   */
  void storeDataInWaveform();

  /// Disable copy construction
  GlobalReadDataContext(const GlobalReadDataContext &copy) = delete;

  /// Disable assignment construction
  GlobalReadDataContext &operator=(const GlobalReadDataContext &assign) = delete;

  /// Move constructor, use the implicitly declared.
  GlobalReadDataContext(GlobalReadDataContext &&) = default;
  GlobalReadDataContext &operator=(GlobalReadDataContext &&) = default;

private:
  static logging::Logger _log;

  /// Waveform wrapped by this GlobalReadDataContext.
  time::PtrWaveform _waveform;
};

} // namespace impl
} // namespace precice
