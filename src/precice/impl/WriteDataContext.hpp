#pragma once

#include <Eigen/Core>

#include "DataContext.hpp"
#include "logging/Logger.hpp"

namespace precice {
namespace impl {

/**
 * @brief Stores one Data object with related mesh. Additionally stores Waveform associated with _providedData.
 *
 * Derived from DataContext
 */
class WriteDataContext : public DataContext {
public:
  WriteDataContext(
      mesh::PtrData data,
      mesh::PtrMesh mesh);

  /**
   * @brief links a write mapping and the mesh context the write mapping requires to this data context
   *
   * @param[in] mappingContext provides context of write mapping
   * @param[in] fromMeshContext provides context of mesh this write mapping is mapping from (_fromData)
   */
  void configureForWriteMapping(MappingContext mappingContext, MeshContext toMeshContext);

private:
  logging::Logger _log{"impl::WriteDataContext"};
};

} // namespace impl
} // namespace precice