#pragma once

#include "mesh/SharedPointer.hpp"
#include "logging/Logger.hpp"
#include "xml/XMLTag.hpp"
#include <vector>
#include <list>
#include <map>
#include <string>

namespace precice {
  namespace mesh {
    class DataConfiguration;
  }
}

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace mesh {

class MeshConfiguration : public xml::XMLTag::Listener
{
public:

  /// Constructor, takes a valid data configuration as argument.
  MeshConfiguration (
    xml::XMLTag&       parent,
    PtrDataConfiguration config );

  void setDimensions ( int dimensions );

  /**
   * @brief Has to be called after parsing all mesh tags.
   *
   * To separate the sub ID setting from the setting of the main mesh IDs makes
   * it possible to have all main IDs in a continously increasing sequence.
   */
  void setMeshSubIDs();

  /// Returns all configured meshes.
  const std::vector<PtrMesh>& meshes() const;

  /// Returns all configured meshes.
  std::vector<PtrMesh>& meshes();

  /// Returns the configured mesh with given name, or NULL.
  mesh::PtrMesh getMesh ( const std::string& meshName ) const;

  virtual void xmlTagCallback(const xml::ConfigurationContext& context, xml::XMLTag& callingTag);

  virtual void xmlEndTagCallback(
          const xml::ConfigurationContext& context,
          xml::XMLTag& callingTag);

  const PtrDataConfiguration& getDataConfiguration() const;

  void addMesh ( const mesh::PtrMesh& mesh );

  std::map<std::string, std::vector<std::string> >& getNeededMeshes(){
    return _neededMeshes;
  }

  void addNeededMesh(
    const std::string& participant,
    const std::string& mesh);

private:

  logging::Logger _log{"mesh::MeshConfiguration"};

  const std::string TAG;
  const std::string ATTR_NAME;
  const std::string ATTR_FLIP_NORMALS;
  const std::string TAG_DATA;
  const std::string TAG_SUB_ID;
  const std::string ATTR_SIDE_INDEX;

  int _dimensions;

  /// Data configuration.
  PtrDataConfiguration _dataConfig;

  /// Configured meshes.
  std::vector<PtrMesh> _meshes;

  bool _setMeshSubIDs;

  std::vector<std::list<std::string> > _meshSubIDs;

  /// to check later if all meshes that any coupling scheme needs are actually used by the participants
  std::map<std::string,std::vector<std::string> > _neededMeshes;
};

}} // namespace precice, mesh

