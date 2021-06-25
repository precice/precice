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

  std::string getParticipantDataName() const;

  int getParticipantDataID() const;

  mesh::PtrData fromData();

  std::string getFromDataName() const;

  int getFromDataID() const;

  mesh::PtrData toData();

  std::string getToDataName() const;

  int getToDataID() const;

  std::string getMeshName() const;

  int getMeshID() const;

  void setMapping(MappingContext mappingContext, mesh::PtrData fromData, mesh::PtrData toData);

  MappingContext _mappingContext;

  bool hasMapping() const;

private:
  mesh::PtrMesh _mesh;

  mesh::PtrData _participantData; // data this participant will write to and read from

  mesh::PtrData _fromData;

  mesh::PtrData _toData;
};

} // namespace impl
} // namespace precice
