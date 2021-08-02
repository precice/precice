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
 *   One of the two must be equal to provdedData.
 * - If this dataContext is not associated with a mapping, fromData and toData will be unset.
 */
class DataContext {

public:
  DataContext(mesh::PtrData data, mesh::PtrMesh mesh);

  mesh::PtrData providedData();

  std::string getDataName() const;

  int getProvidedDataID() const;

  mesh::PtrData fromData();

  int getFromDataID() const;

  mesh::PtrData toData();

  int getToDataID() const;

  std::string getMeshName() const;

  int getMeshID() const;

  void configureForReadMapping(MappingContext mappingContext, MeshContext meshContext);

  void configureForWriteMapping(MappingContext mappingContext, MeshContext meshContext);

  bool hasMapping() const;

  const MappingContext mappingContext() const;

private:
  mesh::PtrMesh _mesh;

  // data this participant will write to and read from
  mesh::PtrData _providedData;

  mesh::PtrData _fromData;

  mesh::PtrData _toData;

  MappingContext _mappingContext;

  bool _hasMapping = false;

  void setMapping(MappingContext mappingContext, mesh::PtrData fromData, mesh::PtrData toData);
};

} // namespace impl
} // namespace precice
