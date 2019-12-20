#pragma once

#include "mapping/SharedPointer.hpp"
#include "mapping/config/MappingConfiguration.hpp"

namespace precice {
namespace impl {

/// Holds a data mapping and related information.
struct MappingContext {
  /// Data mapping.
  mapping::PtrMapping mapping;

  /// id of mesh from which is mapped
  int fromMeshID = -1;

  /// id of mesh to which is mapped
  int toMeshID = -1;

  /// Time of execution of mapping.
  mapping::MappingConfiguration::Timing timing = mapping::MappingConfiguration::INITIAL;

  /// True, if computation and mapping is done repeatedly for single values.
  //bool isIncremental;

  /// True, if data has been mapped already.
  bool hasMappedData = false;
};

} // namespace impl
} // namespace precice
