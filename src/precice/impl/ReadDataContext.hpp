#pragma once

#include <Eigen/Core>

#include "DataContext.hpp"
#include "logging/Logger.hpp"

namespace precice {
namespace impl {

/**
 * @brief Stores one Data object with related mesh.
 *
 * Derived from DataContext
 */
class ReadDataContext : public DataContext {
public:
  ReadDataContext(
      mesh::PtrData data,
      mesh::PtrMesh mesh);

  /**
   * @brief links a read mapping and the mesh context the read mapping requires to this data context
   *
   * @param[in] mappingContext provides context of read mapping
   * @param[in] fromMeshContext provides context of mesh this read mapping is mapping to (_toData)
   */
  void configureForReadMapping(MappingContext mappingContext, MeshContext fromMeshContext);

private:
  logging::Logger _log{"impl::ReadDataContext"};
};

} // namespace impl
} // namespace precice