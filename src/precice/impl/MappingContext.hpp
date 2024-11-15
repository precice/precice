#pragma once

#include "mapping/Mapping.hpp"
#include "mapping/SharedPointer.hpp"
#include "mapping/config/MappingConfiguration.hpp"
namespace precice {
namespace impl {

/// Holds a data mapping and related information.
struct MappingContext {
  /// Data mapping.
  mapping::PtrMapping mapping;

  /// id of mesh from which is mapped
  MeshID fromMeshID = -1;

  /// id of mesh to which is mapped
  MeshID toMeshID = -1;

  /// data which is mapped from mesh
  mesh::PtrData fromData = nullptr;

  /// data which is mapped to mesh
  mesh::PtrData toData = nullptr;

  /// used the automatic rbf alias tag in order to set the mapping
  bool configuredWithAliasTag = false;

  /// Enables gradient data in the corresponding 'from' data class
  void requireGradientData(const std::string &dataName)
  {
    mapping->getInputMesh()->data(dataName)->requireDataGradient();
  }
};

} // namespace impl
} // namespace precice
