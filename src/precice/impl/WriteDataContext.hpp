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
   * @brief Links a MappingContext and the MeshContext required by the write to this WriteDataContext.
   *
   * A write mapping maps _providedData to _toData. A WriteDataContext already has _providedData, but additionally requires _toData.
   * 
   * @param[in] mappingContext provides context of write mapping
   * @param[in] meshContext provides context of mesh this write mapping is mapping to (_toData)
   */
  void configureMapping(MappingContext mappingContext, MeshContext meshContext) override;

private:
  logging::Logger _log{"impl::WriteDataContext"};
};

} // namespace impl
} // namespace precice