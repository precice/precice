#pragma once

#include "MappingContext.hpp"
#include "mesh/SharedPointer.hpp"
#include "mapping/SharedPointer.hpp"

namespace precice {
namespace impl {

/**
 * @brief Stores one Data object with related context. If this dataContext is not associated with a mapping,
 * fromData and toData refer to the same data object.
 */
struct DataContext
{
  bool used;

  mesh::PtrData fromData;

  mesh::PtrData toData;

  mesh::PtrMesh mesh;

  MappingContext mappingContext;

  DataContext():
    used(false),
    fromData(),
    toData(),
    //mesh(),
    mappingContext()
  {}
};

}} // namespace precice, impl
