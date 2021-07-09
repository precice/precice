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

  mesh::PtrData participantData();

  std::string getDataName() const;

  int getParticipantDataID() const;

  mesh::PtrData fromData();

  int getFromDataID() const;

  mesh::PtrData toData();

  int getToDataID() const;

  std::string getMeshName() const;

  int getMeshID() const;

  void setMapping(MappingContext mappingContext, mesh::PtrData fromData, mesh::PtrData toData);

  bool hasMapping() const;

  const MappingContext mappingContext() const;

private:
  mesh::PtrMesh _mesh;

  mesh::PtrData _participantData; // data this participant will write to and read from

  mesh::PtrData _fromData;

  mesh::PtrData _toData;

  MappingContext _mappingContext;

  bool _hasMapping;
};

} // namespace impl
} // namespace precice
