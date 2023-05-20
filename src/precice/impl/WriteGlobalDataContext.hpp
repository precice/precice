#pragma once

#include <Eigen/Core>

#include "DataContext.hpp"
#include "logging/Logger.hpp"
#include "time/Sample.hpp"

namespace precice {
namespace impl {

/**
 * @brief Stores one Global (meshless) Data object. Context stores data to be written to.
 *
 * Derived from DataContext
 */
class WriteGlobalDataContext : public DataContext {
public:
  /**
   * @brief Construct a new WriteGlobalDataContext object without a mapping.
   *
   * @param data Data associated with this WriteGlobalDataContext.
   */
  WriteGlobalDataContext(
      mesh::PtrData data);

  /**
   * @brief Resets provided data, writeDataBuffer, and (optionally) storage.
   *
   * @param atEndOfWindow if true, also the Storage will be reset (useful at end of window to clear storage).
   */
  void resetData(bool atEndOfWindow);

  /**
   * @brief Store value in _writeDataBuffer
   *
   * @param[in] value value of data
   */
  void writeValueIntoDataBuffer(::precice::span<const double> values);

  /**
   * @brief Store data from _writeDataBuffer in persistent storage
   *
   * @param[in] currentTime time data should be associated with
   */
  void storeBufferedData(double currentTime);

  /// Disable copy construction
  WriteGlobalDataContext(const WriteGlobalDataContext &copy) = delete;

  /// Disable assignment construction
  WriteGlobalDataContext &operator=(const WriteGlobalDataContext &assign) = delete;

  /// Move constructor, use the implicitly declared.
  WriteGlobalDataContext(WriteGlobalDataContext &&) = default;
  WriteGlobalDataContext &operator=(WriteGlobalDataContext &&) = default;

private:
  static logging::Logger _log;

  /// @brief Buffer to store written data until it is copied to _providedData->timeStepsStorage()
  time::Sample _writeDataBuffer;
};

} // namespace impl
} // namespace precice
