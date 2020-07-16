#pragma once

#include <list>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#include "logging/Logger.hpp"
#include "mesh/SharedPointer.hpp"
#include "utils/ManageUniqueIDs.hpp"
#include "xml/XMLTag.hpp"

namespace precice {
namespace mesh {
class DataConfiguration;
}
} // namespace precice

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace mesh {

class MeshConfiguration : public xml::XMLTag::Listener {
public:
  /// Constructor, takes a valid data configuration as argument.
  MeshConfiguration(
      xml::XMLTag &        parent,
      PtrDataConfiguration config);

  void setDimensions(int dimensions);

  /// Returns all configured meshes.
  const std::vector<PtrMesh> &meshes() const;

  /// Returns all configured meshes.
  std::vector<PtrMesh> &meshes();

  /// Returns the configured mesh with given name, or NULL.
  mesh::PtrMesh getMesh(const std::string &meshName) const;

  virtual void xmlTagCallback(const xml::ConfigurationContext &context, xml::XMLTag &callingTag);

  virtual void xmlEndTagCallback(
      const xml::ConfigurationContext &context,
      xml::XMLTag &                    callingTag);

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

private:
  logging::Logger _log{"mesh::MeshConfiguration"};

  const std::string TAG;
  const std::string ATTR_NAME;
  const std::string ATTR_FLIP_NORMALS;
  const std::string TAG_DATA;
  const std::string ATTR_SIDE_INDEX;

  int _dimensions;

  /// Data configuration.
  PtrDataConfiguration _dataConfig;

  /// Configured meshes.
  std::vector<PtrMesh> _meshes;

  /// to check later if all meshes that any coupling scheme needs are actually used by the participants
  std::map<std::string, std::vector<std::string>> _neededMeshes;

  std::unique_ptr<utils::ManageUniqueIDs> _meshIdManager;
};

} // namespace mesh
} // namespace precice
