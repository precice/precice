#pragma once

#include <Eigen/Core>

#include "DataContext.hpp"
#include "logging/Logger.hpp"

namespace precice {
namespace impl {

/**
 * @brief Stores one Global (meshless) Data object. Context stores data to be written to.
 *
 * Derived from DataContext
 */
class GlobalWriteDataContext : public DataContext {
public:
  /**
   * @brief Construct a new GlobalWriteDataContext object without a mapping.
   *
   * @param data Data associated with this GlobalWriteDataContext.
   */
  GlobalWriteDataContext(
      mesh::PtrData data);

  /**
   * @brief Get _providedData member.
   *
   * @return mesh::PtrData _providedData.
   */
  mesh::PtrData providedData();

  /// Disable copy construction
  GlobalWriteDataContext(const GlobalWriteDataContext &copy) = delete;

  /// Disable assignment construction
  GlobalWriteDataContext &operator=(const GlobalWriteDataContext &assign) = delete;

  /// Move constructor, use the implicitly declared.
  GlobalWriteDataContext(GlobalWriteDataContext &&) = default;
  GlobalWriteDataContext &operator=(GlobalWriteDataContext &&) = default;

private:
  static logging::Logger _log;
};

} // namespace impl
} // namespace precice
