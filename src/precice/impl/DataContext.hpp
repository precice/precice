#pragma once

#include <map>
#include <optional>
#include <string>

#include "MappingContext.hpp"
#include "MeshContext.hpp"
#include "mapping/SharedPointer.hpp"
#include "mapping/config/MappingConfiguration.hpp"
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
   * @brief Attach a just-in-time mapping to this data context and setup a corresponding MappingDataCache
   *
   * A just-in-time mapping ( \p justInTimeMapping ) might be shared across multiple data contexts (it's a PtrMapping).
   * The MappingDataCache ( \p mappingCache ), however, is unique for each data context.
   *
   * There is no need to put this function into a derived class:
   * The Mapping class knows the direction of the mapping and the DataContext is directly called
   * via the API function. Since all mappings within the DataContexts (stored in a std::vector
   * as \p _mappingContexts ) only call and steer the mapData() execution (and not the
   * computeMapping() execution), we don't register the just-in-time mapping in the \p _mappingContexts
   * data structures. For just-in-time mappings, only the computeMapping() part is relevant, whereas
   * the mapData() part is essentially handled through the API function calls.
   *
   * @param[in] mappingContext The MappingContext which holds the mapping configuration
   * @param[in] meshContext The MeshContext holding the relevant data to map
   */
  void addJustInTimeMapping(MappingContext &mappingContext, MeshContext &meshContext);

  /**
   * @brief Initializes the MappingDataCache
   *
   * Essentially a forward to Mapping::initializeMappingDataCache,
   * where we provide the cache from this class itself.
   *
   * See Mapping::initializeMappingDataCache and the MappingDataCache documentation
   * for more information.
   */
  void initializeMappingDataCache()
  {
    if (mappingCache) {
      justInTimeMapping->initializeMappingDataCache(*mappingCache.get());
      mappingCache->resetData();
    }
  }

  /**
   * @brief Resets the time stamp of the MappingDataCache and potentially resets the data it holds
   *
   * See also the impl::MappingDataCache for more details.
   *
   * @param resetData whether to reset the data or not
   */
  void invalidateMappingCache(bool resetData)
  {
    if (mappingCache) {
      mappingCache->resetTimeStamp();
      if (resetData) {
        mappingCache->resetData();
      }
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
  /// For direct and just-in-time mapping, this is data from the received mesh
  mesh::PtrData _providedData;

  /**
   * @brief Cache for just-in-time mapping.
   *
   * Purpose: For just-in-time mapping, we want to pre-compute as much data as possible to not have
   * expensive operations repeated in every read/write call of the mapping class.
   *
   * The cache itself is updated in the associated Mapping class
   * where it is used. However, storing it in the mapping class is not useful, as we might map
   * multiple data from the same mesh. In preCICE, this leads to one (shared) MappingContext in
   * different DataContexts. Hence, we store the data in this unique Mesh-Data pair.
   *
   * The \p _mappingContexts in this class are a std::vector for multiple mappings associated
   * to the same DataContext. However, for one DataContext (read or write) there can only be one
   * just-in-time mapping.
   * For conventional mappings, multiple MappingContexts can only occur in write direction (write
   * to one local mesh, map it to different remote meshes). Since the just-in-time mapping operates
   * on the remote meshes, this multiplicity cannot occur. Thus, one cache per DataContext is enough.
   * See the documentation of impl::MappingDataCache for more information.
   */
  std::unique_ptr<mapping::impl::MappingDataCache> mappingCache;

  /// The just-in-time mapping for this data context
  mapping::PtrMapping justInTimeMapping;

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
