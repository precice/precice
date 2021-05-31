#pragma once

#include <string>
#include "MappingContext.hpp"
#include "mesh/SharedPointer.hpp"

namespace precice {
namespace impl {

/**
 * @brief Stores one Data object with related context. If this dataContext is not associated with a mapping,
 * fromData and toData refer to the same data object.
 */
struct DataContext {
  bool used = false;

  std::string getName() const;

  mesh::PtrData fromData;

  mesh::PtrData toData;

  mesh::PtrMesh mesh;

  MappingContext mappingContext;
};

} // namespace impl
} // namespace precice
