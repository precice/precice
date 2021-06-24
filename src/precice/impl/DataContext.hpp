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
class DataContext {

public:
  DataContext(mesh::PtrData data, mesh::PtrMesh mesh);

  std::string getFromDataName() const;

  int getFromDataID() const;

  std::string getToDataName() const;

  int getToDataID() const;

  std::string getMeshName() const;

  int getMeshID() const;

  mesh::PtrData fromData();

  void setFromData(mesh::PtrData data);

  mesh::PtrData toData();

  void setToData(mesh::PtrData data);

  MappingContext mappingContext;

private:
  mesh::PtrMesh _mesh;

  mesh::PtrData _fromData;

  mesh::PtrData _toData;
};

} // namespace impl
} // namespace precice
