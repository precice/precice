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
  DataContext(mesh::PtrData data, mesh::PtrMesh mesh);

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

  // links a read mapping to this data context
  void configureForReadMapping(MappingContext mappingContext, MeshContext meshContext);

  // links a write mapping to this data context
  void configureForWriteMapping(MappingContext mappingContext, MeshContext meshContext);

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
