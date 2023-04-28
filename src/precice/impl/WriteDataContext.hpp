#pragma once

#include <Eigen/Core>

#include "DataContext.hpp"
#include "logging/Logger.hpp"

namespace precice {
namespace impl {

/**
 * @brief Stores one Data object with related mesh. Context stores data to be written to and potentially provides a write mapping.
 *
 * Derived from DataContext
 */
class WriteDataContext : public DataContext {
public:
  /**
   * @brief Construct a new WriteDataContext object without a mapping.
   *
   * @param data Data associated with this WriteDataContext.
   * @param mesh Mesh associated with this WriteDataContext.
   */
  WriteDataContext(
      mesh::PtrData data,
      mesh::PtrMesh mesh);

  /**
   * @brief Get _providedData member.
   *
   * @return mesh::PtrData _providedData.
   */
  mesh::PtrData providedData();

  /**
   * @brief Store values in _providedData.values()
   *
   * @param[in] indices ids of data
   * @param[in] values values of data
   */
  void writeValues(const std::vector<int> &indices, const Eigen::Map<const Eigen::VectorXd> values);

  /**
   * @brief Store gradients in _providedData.gradientValues()
   *
   * @param[in] indices ids of data
   * @param[in] gradients gradients of data
   */
  void writeGradientValues(const std::vector<int> &indices, const Eigen::Map<const Eigen::MatrixXd> gradients);

  /**
   * @brief Adds a MappingContext and the MeshContext required by the write mapping to the corresponding WriteDataContext data structures.
   *
   * A write mapping maps _providedData to _toData. A WriteDataContext already has _providedData, but additionally requires _toData.
   *
   * @param[in] mappingContext Context of write mapping
   * @param[in] meshContext Context of mesh this write mapping is mapping to (_toData)
   */
  void appendMappingConfiguration(MappingContext &mappingContext, const MeshContext &meshContext) override;

  /// Disable copy construction
  WriteDataContext(const WriteDataContext &copy) = delete;

  /// Disable assignment construction
  WriteDataContext &operator=(const WriteDataContext &assign) = delete;

  /// Move constructor, use the implicitly declared.
  WriteDataContext(WriteDataContext &&) = default;
  WriteDataContext &operator=(WriteDataContext &&) = default;

private:
  static logging::Logger _log;
};

} // namespace impl
} // namespace precice
