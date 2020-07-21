#include "MeshConfiguration.hpp"
#include <sstream>
#include <stdexcept>
#include "logging/LogMacros.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/config/DataConfiguration.hpp"
#include "utils/Helpers.hpp"
#include "utils/assertion.hpp"
#include "xml/ConfigParser.hpp"
#include "xml/XMLAttribute.hpp"

namespace precice {
namespace mesh {

MeshConfiguration::MeshConfiguration(
    xml::XMLTag &        parent,
    PtrDataConfiguration config)
    : TAG("mesh"),
      ATTR_NAME("name"),
      ATTR_FLIP_NORMALS("flip-normals"),
      TAG_DATA("use-data"),
      ATTR_SIDE_INDEX("side"),
      _dimensions(0),
      _dataConfig(config),
      _meshes(),
      _neededMeshes(),
      _meshIdManager(new utils::ManageUniqueIDs())
{
  using namespace xml;
  std::string doc;
  XMLTag      tag(*this, TAG, xml::XMLTag::OCCUR_ONCE_OR_MORE);
  doc = "Surface mesh consisting of vertices and (optional) of edges and ";
  doc += "triangles (only in 3D). The vertices of a mesh can carry data, ";
  doc += "configured by tag <use-data>. The mesh coordinates have to be ";
  doc += "defined by a participant (see tag <use-mesh>).";
  tag.setDocumentation(doc);

  auto attrName = XMLAttribute<std::string>(ATTR_NAME)
                      .setDocumentation("Unique name for the mesh.");
  tag.addAttribute(attrName);

  auto attrFlipNormals = makeXMLAttribute(ATTR_FLIP_NORMALS, false)
                             .setDocumentation("Flips mesh normal vector directions.");
  tag.addAttribute(attrFlipNormals);

  XMLTag subtagData(*this, TAG_DATA, XMLTag::OCCUR_ARBITRARY);
  doc = "Assigns a before defined data set (see tag <data>) to the mesh.";
  subtagData.setDocumentation(doc);
  attrName.setDocumentation("Name of the data set.");
  subtagData.addAttribute(attrName);
  tag.addSubtag(subtagData);

  parent.addSubtag(tag);
}

void MeshConfiguration::setDimensions(
    int dimensions)
{
  PRECICE_TRACE(dimensions);
  PRECICE_ASSERT((dimensions == 2) || (dimensions == 3), dimensions);
  _dimensions = dimensions;
}

void MeshConfiguration::xmlTagCallback(
    const xml::ConfigurationContext &context,
    xml::XMLTag &                    tag)
{
  PRECICE_TRACE(tag.getName());
  if (tag.getName() == TAG) {
    PRECICE_ASSERT(_dimensions != 0);
    std::string name        = tag.getStringAttributeValue(ATTR_NAME);
    bool        flipNormals = tag.getBooleanAttributeValue(ATTR_FLIP_NORMALS);
    PRECICE_ASSERT(_meshIdManager);
    _meshes.push_back(PtrMesh(new Mesh(name, _dimensions, flipNormals, _meshIdManager->getFreeID())));
  } else if (tag.getName() == TAG_DATA) {
    std::string name  = tag.getStringAttributeValue(ATTR_NAME);
    bool        found = false;
    for (const DataConfiguration::ConfiguredData &data : _dataConfig->data()) {
      if (data.name == name) {
        _meshes.back()->createData(data.name, data.dimensions);
        found = true;
        break;
      }
    }
    if (not found) {
      PRECICE_ERROR("Data with name \"" << name << "\" used by "
                                        << "mesh \"" << _meshes.back()->getName() << "\" is not defined. "
                                        << "Please define a data tag with name=\"" << name << "\".");
    }
  }
}

void MeshConfiguration::xmlEndTagCallback(
    const xml::ConfigurationContext &context,
    xml::XMLTag &                    tag)
{
}

const PtrDataConfiguration &MeshConfiguration::getDataConfiguration() const
{
  return _dataConfig;
}

void MeshConfiguration::addMesh(
    const mesh::PtrMesh &mesh)
{
  for (PtrData dataNewMesh : mesh->data()) {
    bool found = false;
    for (const DataConfiguration::ConfiguredData &data : _dataConfig->data()) {
      if ((dataNewMesh->getName() == data.name) && (dataNewMesh->getDimensions() == data.dimensions)) {
        found = true;
        break;
      }
    }
    PRECICE_ASSERT(found, "Data " << dataNewMesh->getName() << " is not defined. Please define a data tag with name=\"" << dataNewMesh->getName() << "\".");
  }
  _meshes.push_back(mesh);
}

const std::vector<PtrMesh> &MeshConfiguration::meshes() const
{
  return _meshes;
}

std::vector<PtrMesh> &MeshConfiguration::meshes()
{
  return _meshes;
}

mesh::PtrMesh MeshConfiguration::getMesh(
    const std::string &meshName) const
{
  for (const mesh::PtrMesh &mesh : _meshes) {
    if (mesh->getName() == meshName) {
      return mesh;
    }
  }
  return mesh::PtrMesh();
}

void MeshConfiguration::addNeededMesh(
    const std::string &participant,
    const std::string &mesh)
{
  PRECICE_TRACE(participant, mesh);
  if (_neededMeshes.count(participant) == 0) {
    std::vector<std::string> meshes;
    meshes.push_back(mesh);
    _neededMeshes.insert(std::pair<std::string, std::vector<std::string>>(participant, meshes));
  } else if (not utils::contained(mesh, _neededMeshes.find(participant)->second)) {
    _neededMeshes.find(participant)->second.push_back(mesh);
  }
}

} // namespace mesh
} // namespace precice
