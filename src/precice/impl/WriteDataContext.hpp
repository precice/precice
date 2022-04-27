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
   * @brief Links a MappingContext and the MeshContext required by the write to this WriteDataContext.
   *
   * A write mapping maps _providedData to _toData. A WriteDataContext already has _providedData, but additionally requires _toData.
   * 
   * @param[in] mappingContext Context of write mapping
   * @param[in] meshContext Context of mesh this write mapping is mapping to (_toData)
   */
  void addMappingConfiguration(const MappingContext &mappingContext, const MeshContext &meshContext) override;

private:
  static logging::Logger _log;
};

} // namespace impl
} // namespace precice
