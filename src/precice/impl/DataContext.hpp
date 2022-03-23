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
 * - If a DataContext is associated with a mapping, fromData and toData will be set correspondingly.
 *   One of the two must be equal to providedData. fromData and toData must be different.
 * - If a DataContext is not associated with a mapping, fromData and toData will be unset.
 */
class DataContext {
  friend class testing::DataContextFixture; // Make the fixture friend of this class
public:
  /**
   * @brief Get _providedData member.
   *
   * @return mesh::PtrData _providedData.
   */
  mesh::PtrData providedData();

  /**
   * @brief Get the Name of _providedData.
   *
   * @return std::string Name of _providedData.
   */
  std::string getDataName() const;

  /**
   * @brief Get the ID of _providedData.
   *
   * @return int ID of _providedData.
   */
  int getProvidedDataID() const;

  /**
   * @brief Get the ID of _fromData. Used for performing the mapping outside of this class.
   *
   * @return int ID of _fromData.
   */
  int getFromDataID() const;

  /**
   * @brief Purpose unclear. See also https://github.com/precice/precice/issues/1156.
   */
  void resetProvidedData();

  /**
   * @brief Resets _toData to zero. Used before mapping is performed.
   */
  void resetToData();

  /**
   * @brief Get the ID of _toData. Used for performing the mapping outside of this class.
   *
   * @return int ID of _toData.
   */
  int getToDataID() const;

  /**
   * @brief Get the dimensions of _providedData.
   *
   * @return int Dimensions of _providedData.
   */
  int getDataDimensions() const;

  /**
   * @brief Get the name of _mesh.
   *
   * @return std::string Name of _mesh.
   */
  std::string getMeshName() const;

  /**
   * @brief Get the ID of _mesh.
   *
   * @return int ID of _mesh.
   */
  int getMeshID() const;

  /**
   * @brief Informs the user whether this DataContext has a _mappingContext.
   *
   * @return True, if this DataContext is associated with a mapping. False, if not.
   */
  bool hasMapping() const;

  /**
   * @brief Check whether mapping has to be performed.
   *
   * Checks whether a mapping exists for this context and the timing configuration.
   *
   * @return True, if a mapping has to be performed.
   */
  bool isMappingRequired();

  /**
   * @brief Get the _mappingContext associated with this DataContext.
   *
   * @return const MappingContext The _mappingContext of this DataContext.
   */
  const MappingContext mappingContext() const;

  /**
   * @brief Informs the user whether this DataContext has a read mapping.
   *
   * @return True, if DataContext has a read mapping.
   */
  bool hasReadMapping() const;

  /**
   * @brief Informs the user whether this DataContext has a write mapping.
   *
   * @return True, if DataContext has a write mapping.
   */
  bool hasWriteMapping() const;

  /**
   * @brief Links a MappingContext and the MeshContext required by the mapping to this DataContext.
   *
   * A mapping maps the given data from or to _providedData (depending on whether it is a read or write mapping).
   *
   * @param[in] mappingContext Context of the mapping
   * @param[in] meshContext Context of mesh this mapping is mapping from or to
   */
  virtual void configureMapping(const MappingContext &mappingContext, const MeshContext &meshContext) = 0;

protected:
  /**
   * @brief Construct a new DataContext without a mapping. Protected, because only ReadDataContext and WriteDataContext should use this constructor.
   *
   * @param data Data associated with this DataContext.
   * @param mesh Mesh associated with this DataContext.
   */
  DataContext(mesh::PtrData data, mesh::PtrMesh mesh);

  /// Defines the mapping associated to this DataContext. A DataContext may also exist without a mapping.
  MappingContext _mappingContext;

  /// Data this participant will write to and read from
  mesh::PtrData _providedData;

  /// If a mapping exists, mesh::PtrData the mapping maps from.
  mesh::PtrData _fromData;

  /// If a mapping exists, mesh::PtrData the mapping maps from.
  mesh::PtrData _toData;

  /**
   * @brief Helper to set _mappingContext, _fromData and _toData.
   *
   * @param mappingContext MappingContext this DataContext will be associated to.
   * @param fromData Data the mapping maps from.
   * @param toData Data the mapping maps to.
   */
  void setMapping(MappingContext mappingContext, mesh::PtrData fromData, mesh::PtrData toData);

private:
  /// Mesh associated with _providedData.
  mesh::PtrMesh _mesh;

  static logging::Logger _log;
};

} // namespace impl
} // namespace precice
