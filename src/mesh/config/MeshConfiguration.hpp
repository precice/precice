#pragma once

#include <list>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#include "logging/Logger.hpp"
#include "mesh/Data.hpp"
#include "mesh/SharedPointer.hpp"
#include "utils/ManageUniqueIDs.hpp"
#include "xml/XMLTag.hpp"

namespace precice::mesh {
class DataConfiguration;
}

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice::mesh {

class MeshConfiguration : public xml::XMLTag::Listener {
public:
  /// Constructor, takes a valid data configuration as argument.
  MeshConfiguration(
      xml::XMLTag &        parent,
      PtrDataConfiguration config);

  /// Returns all configured meshes.
  const std::vector<PtrMesh> &meshes() const;

  /// Returns all configured meshes.
  std::vector<PtrMesh> &meshes();

  /// Returns whether Mesh has Data with the dataName
  bool hasMeshName(const std::string &meshName) const;

  /// Returns the configured mesh with given name, or NULL.
  mesh::PtrMesh getMesh(const std::string &meshName) const;

  static mesh::PtrMesh getIndirectAccessMesh(int dimension);

  void xmlTagCallback(const xml::ConfigurationContext &context, xml::XMLTag &callingTag) override;

  void xmlEndTagCallback(
      const xml::ConfigurationContext &context,
      xml::XMLTag &                    callingTag) override;

  const PtrDataConfiguration &getDataConfiguration() const;

  void addMesh(const mesh::PtrMesh &mesh);

  std::map<std::string, std::vector<std::string>> &getNeededMeshes()
  {
    return _neededMeshes;
  }

  void addNeededMesh(
      const std::string &participant,
      const std::string &mesh);

  std::unique_ptr<utils::ManageUniqueIDs> extractMeshIdManager()
  {
    return std::move(_meshIdManager);
  }

  /// Initialize the map between meshes and dimensions, for unit tests that directly create mesh objects without going through the config reading.
  void insertMeshToMeshDimensionsMap(const std::string &mesh,
                                     int                dimensions);

private:
  logging::Logger _log{"mesh::MeshConfiguration"};

  const std::string TAG;
  const std::string ATTR_NAME;
  const std::string ATTR_DIMENSIONS;
  const std::string TAG_DATA;
  const std::string ATTR_SIDE_INDEX;

  std::map<std::string, int> _meshDimensionsMap;

  /// Data configuration.
  PtrDataConfiguration _dataConfig;

  /// Get the number of dimensions that data values of this type (scalar/vector) have on this mesh
  int getDataDimensions(const std::string &meshName, const Data::typeName typeName) const;

  /// Configured meshes.
  std::vector<PtrMesh> _meshes;

  /// to check later if all meshes that any coupling scheme needs are actually used by the participants
  std::map<std::string, std::vector<std::string>> _neededMeshes;

  std::unique_ptr<utils::ManageUniqueIDs> _meshIdManager;

  utils::ManageUniqueIDs _dataIDManager;
};

} // namespace precice::mesh
