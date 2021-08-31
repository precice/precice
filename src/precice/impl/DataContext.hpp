#pragma once

#include <string>
#include "MappingContext.hpp"
#include "MeshContext.hpp"
#include "mesh/SharedPointer.hpp"

namespace precice {
namespace impl {

/**
 * @brief Stores one Data object with related mesh.
 *
 * - If this dataContext is associated with a mapping, fromData and toData will be set correspondingly.
 *   One of the two must be equal to providedData. fromData and toData must be different.
 * - If this dataContext is not associated with a mapping, fromData and toData will be unset.
 */
class DataContext {

public:
  DataContext(mesh::PtrData data, mesh::PtrMesh mesh); // @ todo use mesh context for constructor? Would be more consistent with configureForReadMapping and configureForWriteMapping

  mesh::PtrData providedData();

  mesh::PtrData toData();

  std::string getDataName() const;

  int getProvidedDataID() const;

  int getFromDataID() const;

  void resetProvidedData();

  void resetToData();

  int getToDataID() const;

  std::string getMeshName() const;

  int getMeshID() const;

  /**
   * @brief links a read mapping and the mesh context the read mapping requires to this data context
   *
   * @param[in] mappingContext provides context of read mapping
   * @param[in] fromMeshContext provides context of mesh this read mapping is mapping to (_toData)
   */
  void configureForReadMapping(MappingContext mappingContext, MeshContext fromMeshContext);

  /**
   * @brief links a write mapping and the mesh context the write mapping requires to this data context
   *
   * @param[in] mappingContext provides context of write mapping
   * @param[in] fromMeshContext provides context of mesh this write mapping is mapping from (_fromData)
   */
  void configureForWriteMapping(MappingContext mappingContext, MeshContext toMeshContext);

  bool hasMapping() const;

  bool hasReadMapping() const;

  bool hasWriteMapping() const;

  const MappingContext mappingContext() const;

private:
  mesh::PtrMesh _mesh;

  // data this participant will write to and read from
  mesh::PtrData _providedData;

  mesh::PtrData _fromData;

  mesh::PtrData _toData;

  MappingContext _mappingContext;

  // helper function for creating read and write mappings
  void setMapping(MappingContext mappingContext, mesh::PtrData fromData, mesh::PtrData toData);
};

} // namespace impl
} // namespace precice
