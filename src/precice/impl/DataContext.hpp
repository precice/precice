#pragma once

#include <string>
#include "MappingContext.hpp"
#include "MeshContext.hpp"
#include "mesh/SharedPointer.hpp"

namespace precice {

namespace testing {
// Forward declaration to friend the boost test struct
class DataContextFixture;
} // namespace testing

namespace impl {

/**
 * @brief Stores one Data object with related mesh.
 *
 * - If this dataContext is associated with a mapping, fromData and toData will be set correspondingly.
 *   One of the two must be equal to providedData. fromData and toData must be different.
 * - If this dataContext is not associated with a mapping, fromData and toData will be unset.
 */
class DataContext {
  friend class testing::DataContextFixture; // Make the fixture friend of this class
public:
  mesh::PtrData providedData();

  mesh::PtrData toData();

  std::string getDataName() const;

  int getProvidedDataID() const;

  int getFromDataID() const;

  void resetProvidedData();

  void resetToData();

  int getToDataID() const;

  int getDataDimensions() const;

  std::string getMeshName() const;

  int getMeshID() const;

  bool hasMapping() const;

  const MappingContext mappingContext() const;

  bool hasReadMapping() const;

  bool hasWriteMapping() const;

protected:
  DataContext(mesh::PtrData data, mesh::PtrMesh mesh);

  // data this participant will write to and read from
  mesh::PtrData _providedData;

  mesh::PtrData _fromData;

  mesh::PtrData _toData;

  MappingContext _mappingContext;

  /// helper function for creating read and write mappings
  void setMapping(MappingContext mappingContext, mesh::PtrData fromData, mesh::PtrData toData);

private:
  mutable logging::Logger _log{"impl::DataContext"};

  mesh::PtrMesh _mesh;
};

} // namespace impl
} // namespace precice
