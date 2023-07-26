#include "MeshConfiguration.hpp"
#include <memory>

#include <sstream>
#include <stdexcept>
#include <utility>

#include "logging/LogMacros.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/config/DataConfiguration.hpp"
#include "utils/Helpers.hpp"
#include "utils/assertion.hpp"
#include "xml/ConfigParser.hpp"
#include "xml/XMLAttribute.hpp"

namespace precice::mesh {

MeshConfiguration::MeshConfiguration(
    xml::XMLTag &        parent,
    PtrDataConfiguration config)
    : TAG("mesh"),
      ATTR_NAME("name"),
      ATTR_DIMENSIONS("dimensions"),
      TAG_DATA("use-data"),
      ATTR_SIDE_INDEX("side"),
      _meshDimensionsMap(),
      _dataConfig(std::move(config)),
      _meshes(),
      _neededMeshes(),
      _meshIdManager(new utils::ManageUniqueIDs())
{
  using namespace xml;
  std::string doc;
  XMLTag      tag(*this, TAG, xml::XMLTag::OCCUR_ONCE_OR_MORE);
  doc = "Surface mesh consisting of vertices and optional connectivity information. "
        "The vertices of a mesh can carry data, "
        "configured by tags <use-data>. The mesh coordinates have to be "
        "defined by a participant (see tag <provide-mesh>).";
  tag.setDocumentation(doc);

  auto attrName = XMLAttribute<std::string>(ATTR_NAME)
                      .setDocumentation("Unique name for the mesh.");
  tag.addAttribute(attrName);

  auto attrDimensions = XMLAttribute<int>(ATTR_DIMENSIONS)
                            .setDocumentation("Spatial dimensions of mesh")
                            .setOptions({2, 3});
  tag.addAttribute(attrDimensions);

  XMLTag subtagData(*this, TAG_DATA, XMLTag::OCCUR_ARBITRARY);
  doc = "Assigns a before defined data set (see tag <data>) to the mesh.";
  subtagData.setDocumentation(doc);
  attrName.setDocumentation("Name of the data set.");
  subtagData.addAttribute(attrName);
  tag.addSubtag(subtagData);

  parent.addSubtag(tag);
}

void MeshConfiguration::xmlTagCallback(
    const xml::ConfigurationContext &context,
    xml::XMLTag &                    tag)
{
  PRECICE_TRACE(tag.getName());
  if (tag.getName() == TAG) {
    std::string name       = tag.getStringAttributeValue(ATTR_NAME);
    int         dimensions = tag.getIntAttributeValue(ATTR_DIMENSIONS);
    _meshDimensionsMap.insert<std::pair<std::string, int>>(name, dimensions);
    PRECICE_ASSERT(dimensions != 0);
    PRECICE_ASSERT(_meshIdManager);
    _meshes.push_back(std::make_shared<Mesh>(name, dimensions, _meshIdManager->getFreeID()));
  } else if (tag.getName() == TAG_DATA) {
    std::string name  = tag.getStringAttributeValue(ATTR_NAME);
    bool        found = false;
    for (const DataConfiguration::ConfiguredData &data : _dataConfig->data()) {
      auto dataDimensions = getDataDimensions(_meshes.back()->getName(), data.typeName);
      if (data.name == name) {
        _meshes.back()->createData(data.name, dataDimensions, _dataIDManager.getFreeID(), data.waveformDegree);
        found = true;
        break;
      }
    }
    if (not found) {
      PRECICE_ERROR("Data with name \"{}\" used by mesh \"{}\" is not defined. "
                    "Please define a data tag with name=\"{}\".",
                    name, _meshes.back()->getName(), name);
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
  for (const PtrData &dataNewMesh : mesh->data()) {
    bool found = false;
    for (const DataConfiguration::ConfiguredData &data : _dataConfig->data()) {
      if ((dataNewMesh->getName() == data.name) && (dataNewMesh->getDimensions() == getDataDimensions(data.typeName))) {
        found = true;
        break;
      }
    }
    PRECICE_CHECK(found, "Data {0} is not defined. Please define a data tag with name=\"{0}\".", dataNewMesh->getName());
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

bool MeshConfiguration::hasMeshName(const std::string &meshName) const
{
  auto iter = std::find_if(_meshes.begin(), _meshes.end(), [&meshName](const auto &mptr) {
    return mptr->getName() == meshName;
  });
  return iter != _meshes.end(); // if name was not found in _meshes, iter == _meshes.end()
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

int MeshConfiguration::getDataDimensions(const std::string &meshName, const std::string &typeName)
{
  // TODO: Get the values from elsewhere / make an enum
  if (typeName == "vector") {
    return _meshDimensionsMap[meshName];
  } else if (typeName == "scalar") {
    return 1;
  }
  // We should never reach this point
  PRECICE_UNREACHABLE("Unknown data type \"{}\".", typeName);
};

} // namespace precice::mesh
