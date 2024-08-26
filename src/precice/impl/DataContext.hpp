#pragma once

#include <map>
#include <optional>
#include <string>

#include "MappingContext.hpp"
#include "MeshContext.hpp"
#include "mapping/NearestNeighborMapping.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Utils.hpp"
namespace precice {

namespace testing {
// Forward declaration to friend the boost test struct
class DataContextFixture;
} // namespace testing

namespace impl {

/**
 * @brief Stores one Data object with related mesh.
 *
 * - For each mapping that is added to the data context fromData and toData will be set correspondingly.
 *   Either fromData or toData must be equal to providedData. fromData and toData must be different.
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
   * @brief Resets initial guesses of transient mappings to zero.
   */
  void resetInitialGuesses();

  /**
   * @brief Get the dimensions of _providedData.
   *
   * @return int Dimensions of _providedData.
   */
  int getDataDimensions() const;

  /**
   * @brief Get the spatial dimensions of _providedData.
   *
   * @return int Spatial dimensions of _providedData.
   */
  int getSpatialDimensions() const;

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
   * @brief Returns whether _providedData has gradient
   *
   * @return true, if it has gradient
   * @return false, if it has gradient
   */
  bool hasGradient() const;

  /**
   * @brief Perform the mapping for mapping contexts and the corresponding data context (from and to data)
   *
   * @param[in] after only map samples after this optional time
   * @param[in] skipZero set output sample to zero if the input sample is zero too
   *
   * @return the number of performed mappings
   */
  int mapData(std::optional<double> after = std::nullopt, bool skipZero = false);

  /**
   * @brief Adds a MappingContext and the MeshContext required by the mapping to the corresponding DataContext data structures.
   *
   * A mapping maps the given data from or to _providedData (depending on whether it is a read or write mapping).
   *
   * @param[in] mappingContext Context of the mapping
   * @param[in] meshContext Context of mesh this mapping is mapping from or to
   */
  virtual void appendMappingConfiguration(MappingContext &mappingContext, const MeshContext &meshContext) = 0;

  /**
   * @brief Informs the user whether this DataContext has any _mappingContext.
   *
   * @return True, if this DataContext is associated with a mapping. False, if not.
   */
  bool hasMapping() const;

  template <typename Container>
  std::optional<std::size_t> locateInvalidVertexID(const Container &c)
  {
    return mesh::locateInvalidVertexID(*_mesh, c);
  }

  /**
   * @brief
   *
   * No need to put this function into a derived class.
   * The Mapping class knows the direction and the DataContext states read or write.
   * Mappings in the DataContext only affect the mapData() steering.
   * Thus, we don't add the indirect mapping anywhere in the conventional mapping
   * data structures, i.e., _mappingContexts
   *
   * @param mappingContext
   */
  void addIndirectAccessMapping(MappingContext mappingContext, MeshContext meshContext);

  void invalidateMappingCache()
  {
    if (mappingCache) {
      mappingCache->setTimeStamp(-1);
    }
  }

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

  /**
   * @brief Cache for indirect mesh access.
   *
   * Purpose: For indirect access, we want to pre-compute as much data as possible to not have
   * expensive operations repeated in every read/write call of the mapping class.
   *
   * The cache itself is updated in the associated Mapping class
   * where it is used. However, storing it in the mapping class is not useful, as we might map
   * multiple data from the same mesh. In preCICE, this leads to one (shared) MappingContext in
   * different DataContexts. Hence, we store the data in this unique Mesh-Data pair.
   *
   * The \ref _mappingContexts in this class are a std::vector for multiple mappings associated
   * to the same DataContext. However, for one DataContext (read or write) there can only be one
   * indirect mapping/indirect mesh access.
   * For conventional mappings, multiple MappingContexts can only occur in write direction (write
   * to one local mesh, map it to different remote meshes). Since the indirect mapping operates
   * on the remote meshes, this multiplicity cannot occur. Thus, one cache per DataContext is enough.
   */
  std::unique_ptr<mapping::MappingDataCache>       mappingCache;
  std::shared_ptr<mapping::Mapping> indirectMapping;
  /**
   * @brief Helper to append a mappingContext, fromData and toData to the corresponding data containers
   *
   * @param mappingContext MappingContext this DataContext will be associated to.
   *
   * @note Only unique mappings may be appended. In case the same mapping is appended twice, an error is raised.
   */
  void appendMapping(MappingContext mappingContext);

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

  /**
   * @brief Get the number of vertices of mesh
   *
   * @return int number of vertices
   */
  int getMeshVertexCount() const;

  /// Returns true if the given vertexID is valid
  bool isValidVertexID(const VertexID id) const;

  mesh::PtrMesh _mesh;

private:
  /// Unique mesh associated with _providedData.

  static logging::Logger _log;

  using FromToDataIDs = std::pair<int, int>;
  std::map<FromToDataIDs, Eigen::VectorXd> _initialGuesses;
};

} // namespace impl
} // namespace precice
