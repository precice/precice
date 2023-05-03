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
   * @brief Get _providedData member.
   *
   * @return mesh::PtrData _providedData.
   */
  mesh::PtrData providedData();

  /// Disable copy construction
  WriteGlobalDataContext(const WriteGlobalDataContext &copy) = delete;

  /// Disable assignment construction
  WriteGlobalDataContext &operator=(const WriteGlobalDataContext &assign) = delete;

  /// Move constructor, use the implicitly declared.
  WriteGlobalDataContext(WriteGlobalDataContext &&) = default;
  WriteGlobalDataContext &operator=(WriteGlobalDataContext &&) = default;

private:
  static logging::Logger _log;
};

} // namespace impl
} // namespace precice
