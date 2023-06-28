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
class ReadGlobalDataContext : public DataContext {
public:
  /**
   * @brief Construct a new ReadGlobalDataContext object without a mapping.
   *
   * @param data Data associated with this ReadGlobalDataContext.
   */
  ReadGlobalDataContext(
      mesh::PtrData data);

  /**
   * @brief Gets _interpolationOrder of _waveform
   *
   * @return _interpolationOrder of _waveform
   */
  int getInterpolationOrder() const;

  /**
   * @brief Samples data at a given point in time within the current time window and writes it to the given span
   *
   * @param[in] normalizedDt Point in time where waveform is sampled. Must be normalized to [0,1], where 0 refers to the beginning and 1 to the end of the current time window.
   * @param[in] value read data at time normalizedDt will be returned into this span
   */
  void readValue(double normalizedDt, ::precice::span<double> value) const;

  /// Disable copy construction
  ReadGlobalDataContext(const ReadGlobalDataContext &copy) = delete;

  /// Disable assignment construction
  ReadGlobalDataContext &operator=(const ReadGlobalDataContext &assign) = delete;

  /// Move constructor, use the implicitly declared.
  ReadGlobalDataContext(ReadGlobalDataContext &&) = default;
  ReadGlobalDataContext &operator=(ReadGlobalDataContext &&) = default;

private:
  static logging::Logger _log;

  /// Waveform wrapped by this ReadGlobalDataContext.
  time::PtrWaveform _waveform;
};

} // namespace impl
} // namespace precice
