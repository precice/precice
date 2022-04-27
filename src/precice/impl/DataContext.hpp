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
 * - A DataContext can be associated with multiple mappings, fromData and toData
 */
class DataContext {
  friend class testing::DataContextFixture; // Make the fixture friend of this class
public:
  /**
   * @brief Get the Name of _providedData.
   *
   * @return std::string Name of _providedData.
   */
  std::string getDataName() const;

  /**
   * @brief Resets provided data and (if mapping exists) fromData or toData.
   */
  void resetData();

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
  MeshID getMeshID() const;

  /**
   * @brief Check whether any mapping has to be performed.
   *
   * Checks whether any mapping exists for this context and the corresponding timing configuration.
   *
   * @return True, if any mapping has to be performed.
   */
  bool isMappingRequired();

  /**
   * @brief Perform the mapping for all mapping contexts and the corresponding data context (from and to data)
   */
  void mapData();

  /**
   * @brief Adds a MappingContext and the MeshContext required by the mapping to the correspnding DataContext data structures.
   *
   * A mapping maps the given data from or to _providedData (depending on whether it is a read or write mapping).
   *
   * @param[in] mappingContext Context of the mapping
   * @param[in] meshContext Context of mesh this mapping is mapping from or to
   */
  virtual void addMappingConfiguration(const MappingContext &mappingContext, const MeshContext &meshContext) = 0;

protected:
  /**
   * @brief Construct a new DataContext without a mapping. Protected, because only ReadDataContext and WriteDataContext should use this constructor.
   *
   * @param data Data associated with this DataContext.
   * @param mesh Mesh associated with this DataContext.
   */
  DataContext(mesh::PtrData data, mesh::PtrMesh mesh);

  /// Defines all mappings associated to this DataContext. A DataContext may also exist without a mapping.
  std::vector<MappingContext> _mappingContexts;

  /// Unique data this context is associated with
  mesh::PtrData _providedData;

  /// If a mapping exists, collection of mesh::PtrData the mapping maps from.
  std::vector<mesh::PtrData> _fromDatas;

  /// If a mapping exists, collection of mesh::PtrData the mapping maps from.
  std::vector<mesh::PtrData> _toDatas;

  /**
   * @brief Helper to append a _mappingContext, _fromData and _toData to the corresponding data containers
   *
   * @param mappingContext MappingContext this DataContext will be associated to.
   * @param fromData Data the mapping maps from.
   * @param toData Data the mapping maps to.
   *
   * @note Only unique mappings may be appended. In case the same mapping is appended twice, an error is raised.
   */
  void appendMapping(MappingContext mappingContext, mesh::PtrData fromData, mesh::PtrData toData);

  /**
   * @brief Informs the user whether this DataContext has any read mapping.
   *
   * @return True, if DataContext has any read mapping.
   */
  bool hasReadMapping() const;

  /**
   * @brief Informs the user whether this DataContext has any write mapping.
   *
   * @return True, if DataContext has any write mapping.
   */
  bool hasWriteMapping() const;

private:
  /// Unique mesh associated with _providedData.
  mesh::PtrMesh _mesh;

  static logging::Logger _log;

  /**
   * @brief Get the ID of the data in the _fromDatas container. Used for performing the mapping outside of this class.
   *
   * @param[in] dataVectorIndex Index of the '_fromDatas' container this data context holds
   *
   * @return DataID ID of _fromDatas.
   */
  DataID getFromDataID(int dataVectorIndex) const;

  /**
   * @brief Get the ID of the data in the _toDatas container. Used for performing the mapping outside of this class.
   *
   * @param[in] dataVectorIndex Index of the '_toDatas' vector this data context holds
   *
   * @return DataID ID of _toDatas.
   */
  DataID getToDataID(int dataVectorIndex) const;

  /**
   * @brief Informs the user whether this DataContext has any _mappingContext.
   *
   * @return True, if this DataContext is associated with a mapping. False, if not.
   */
  bool hasMapping() const;
};

} // namespace impl
} // namespace precice
