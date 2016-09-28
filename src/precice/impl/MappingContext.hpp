#pragma once

#include "mesh/SharedPointer.hpp"
#include "mapping/SharedPointer.hpp"
#include "mapping/config/MappingConfiguration.hpp"

namespace precice {
namespace impl {

/**
 * @brief Holds a data mapping and related information.
 */
struct MappingContext
{
  // @brief Data mapping.
  mapping::PtrMapping mapping;

  // @brief id of mesh from which is mapped
  int fromMeshID;

  // @brief id of mesh to which is mapped
  int toMeshID;

  // @brief Time of execution of mapping.
  mapping::MappingConfiguration::Timing timing;

  // @brief True, if computation and mapping is done repeatedly for single values.
  //bool isIncremental;

  // @brief True, if data has been mapped already.
  bool hasMappedData;

  /**
   * @brief Constructor.
   */
  MappingContext()
  : mapping(),
    fromMeshID(-1),
    toMeshID(-1),
    timing(mapping::MappingConfiguration::INITIAL),
    //isIncremental(false),
    hasMappedData(false)
  {}
};

}} // namespace precice, impl
