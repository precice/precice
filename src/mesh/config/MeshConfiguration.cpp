// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "MeshConfiguration.hpp"
#include "mesh/config/DataConfiguration.hpp"
#include "mesh/Mesh.hpp"
#include "utils/xml/XMLAttribute.hpp"
#include "utils/Globals.hpp"
#include <sstream>

namespace precice {
namespace mesh {

tarch::logging::Log MeshConfiguration:: _log ( "precice::mesh::MeshConfiguration" );

//const std::string& MeshConfiguration:: getTag ()
//{
//  static std::string tag ( "mesh" );
//  return tag;
//}

MeshConfiguration:: MeshConfiguration
(
  utils::XMLTag&       parent,
  PtrDataConfiguration config )
:
  TAG("mesh"),
  ATTR_NAME("name"),
  ATTR_FLIP_NORMALS("flip-normals"),
  TAG_DATA("use-data"),
  TAG_SPACETREE("use-spacetree"),
  TAG_SUB_ID("sub-id"),
  ATTR_SIDE_INDEX("side"),
  _dimensions(0),
  //_isValid(false),
  _dataConfig(config),
  _meshes(),
  _setMeshSubIDs(),
  _meshSubIDs(),
  _spacetreeNames(),
  _neededMeshes()
{
  using namespace utils;
  std::string doc;
  XMLTag tag(*this, TAG, utils::XMLTag::OCCUR_ONCE_OR_MORE);
  doc = "Surface mesh consisting of vertices and (optional) of edges and ";
  doc += "triangles (only in 3D). The vertices of a mesh can carry data, ";
  doc += "configured by tag <use-data>. The geometry of a mesh is either defined ";
  doc += "by a tag <geometry>, or by a participant (see tag <use-mesh>).";
  tag.setDocumentation(doc);

  XMLAttribute<std::string> attrName(ATTR_NAME);
  attrName.setDocumentation("Unique name for the mesh.");
  tag.addAttribute(attrName);

  XMLAttribute<bool> attrFlipNormals(ATTR_FLIP_NORMALS);
  attrFlipNormals.setDocumentation("Flips mesh normal vector directions.");
  attrFlipNormals.setDefaultValue(false);
  tag.addAttribute(attrFlipNormals);

  XMLTag subtagData(*this, TAG_DATA, XMLTag::OCCUR_ARBITRARY);
  doc = "Assigns a before defined data set (see tag <data>) to the mesh.";
  subtagData.setDocumentation(doc);
  attrName.setDocumentation("Name of the data set.");
  subtagData.addAttribute(attrName);
  tag.addSubtag(subtagData);

  utils::XMLTag tagSubID(*this, TAG_SUB_ID, utils::XMLTag::OCCUR_ARBITRARY);
  doc = "Every mesh has a global ID (determined by preCICE). It is possible ";
  doc += "to set additional sub-ids. This sub-ids are used by certain geomtries ";
  doc += "(see tag <geometry>) to distinguish parts of the mesh in geometry ";
  doc += "queries.";
  tagSubID.setDocumentation(doc);
  utils::XMLAttribute<int> attrSideIndex(ATTR_SIDE_INDEX);
  tagSubID.addAttribute(attrSideIndex);
  tag.addSubtag(tagSubID);

  XMLTag subtagUseSpacetree(*this, TAG_SPACETREE, XMLTag::OCCUR_NOT_OR_ONCE);
  doc = "A defined spacetree (see tag <spacetree>) can be used to accelerate ";
  doc += "spatial queries on the mesh.";
  subtagUseSpacetree.setDocumentation(doc);
  attrName.setDocumentation("Name of the spacetree");
  subtagUseSpacetree.addAttribute(attrName);
  tag.addSubtag(subtagUseSpacetree);

  parent.addSubtag(tag);
}

void MeshConfiguration:: setDimensions
(
  int dimensions )
{
  preciceTrace1("setDimensions()", dimensions);
  assertion1((dimensions == 2) || (dimensions == 3), dimensions);
  _dimensions = dimensions;
}

//bool MeshConfiguration:: parseSubtag
//(
//  utils::XMLTag::XMLReader* xmlReader )
//{
//  preciceTrace("parseSubtag()");
//
//  assertion(not _setMeshSubIDs);
//
//  XMLTag tag ( TAG, XMLTag::OCCUR_ONCE );
//  XMLAttribute<std::string> attrName ( ATTR_NAME );
//  tag.addAttribute ( attrName );
//
//  XMLAttribute<bool> attrFlipNormals ( ATTR_FLIP_NORMALS );
//  attrFlipNormals.setDefaultValue ( false );
//  tag.addAttribute ( attrFlipNormals );
//
//  XMLTag subtagData ( TAG_DATA, XMLTag::OCCUR_ARBITRARY );
//  subtagData.addAttribute ( attrName );
//  tag.addSubtag ( subtagData );
//
//  utils::XMLTag tagSubID ( TAG_SUB_ID, utils::XMLTag::OCCUR_ARBITRARY );
//  utils::XMLAttribute<int> attrSideIndex ( ATTR_SIDE_INDEX );
//  tagSubID.addAttribute ( attrSideIndex );
//  tag.addSubtag ( tagSubID );
//
//  XMLTag subtagUseSpacetree ( TAG_SPACETREE, XMLTag::OCCUR_NOT_OR_ONCE );
//  subtagUseSpacetree.addAttribute ( attrName );
//  tag.addSubtag ( subtagUseSpacetree );
//
//  _isValid = _tag.parse(xmlReader);
//  return _isValid;
//}

void MeshConfiguration:: xmlTagCallback
(
  utils::XMLTag& tag )
{
  preciceTrace1("xmlTagCallback()", tag.getName());
  if (tag.getName() == TAG){
    assertion(_dimensions != 0);
    std::string name = tag.getStringAttributeValue(ATTR_NAME);
    bool flipNormals = tag.getBooleanAttributeValue(ATTR_FLIP_NORMALS);
    _meshes.push_back(PtrMesh(new Mesh(name, _dimensions, flipNormals)));
    _meshSubIDs.push_back(std::list<std::string>());
  }
  else if (tag.getName() == TAG_SUB_ID){
    int side = tag.getIntAttributeValue(ATTR_SIDE_INDEX);
    std::stringstream conv;
    conv << "side-" << side;
    _meshSubIDs.back().push_back(conv.str());
  }
  else if (tag.getName() == TAG_DATA){
    std::string name = tag.getStringAttributeValue(ATTR_NAME);
    bool found = false;
    foreach (const DataConfiguration::ConfiguredData& data, _dataConfig->data()){
      if (data.name == name){
        _meshes.back()->createData(data.name, data.dimensions);
        found = true;
        break;
      }
    }
    if (not found){
      std::ostringstream stream;
      stream << "Data with name \"" << name << "\" is not available during "
             << "configuration of mesh \"" << _meshes.back()->getName() << "\"";
      throw stream.str();
    }
  }
  else if (tag.getName() == TAG_SPACETREE){
    std::string name = tag.getStringAttributeValue(ATTR_NAME);
    std::string meshName ( _meshes.back()->getName() );
    if (utils::contained(meshName, _spacetreeNames)){
      std::ostringstream stream;
      stream << "Mesh \"" << meshName << "\" can only use one spacetree";
      throw stream.str();
    }
    _spacetreeNames[meshName] = name;
  }
}

void MeshConfiguration:: xmlEndTagCallback
(
  utils::XMLTag& tag )
{
}

const PtrDataConfiguration& MeshConfiguration:: getDataConfiguration() const
{
  return _dataConfig;
}

void MeshConfiguration:: addMesh
(
  const mesh::PtrMesh& mesh  )
{
  foreach (PtrData dataNewMesh, mesh->data()){
    bool found = false;
    foreach (const DataConfiguration::ConfiguredData & data, _dataConfig->data()){
      if ((dataNewMesh->getName() == data.name)
          && (dataNewMesh->getDimensions() == data.dimensions))
      {
        found = true;
        break;
      }
    }
    preciceCheck(found, "addMesh()", "Data " << dataNewMesh->getName()
                 << " is not available in data configuration!");
  }
  _meshes.push_back(mesh);
}

void MeshConfiguration:: setMeshSubIDs()
{
  //assertion ( _isValid );
  assertion ( _meshes.size() == _meshSubIDs.size() );
  assertion ( not _setMeshSubIDs );
  for ( size_t i=0; i < _meshes.size(); i++ ) {
    foreach ( const std::string & subIDName, _meshSubIDs[i] ) {
      _meshes[i]->setSubID ( subIDName );
    }
  }
  _setMeshSubIDs = true;
}

//bool MeshConfiguration:: isValid() const
//{
//  return _isValid;
//}

const std::vector<PtrMesh>& MeshConfiguration:: meshes() const
{
  return _meshes;
}

std::vector<PtrMesh>& MeshConfiguration:: meshes()
{
  return _meshes;
}

mesh::PtrMesh MeshConfiguration:: getMesh
(
  const std::string& meshName ) const
{
  foreach ( const mesh::PtrMesh & mesh, _meshes ) {
    if ( mesh->getName() == meshName ) {
      return mesh;
    }
  }
  return mesh::PtrMesh();
}

bool MeshConfiguration:: doesMeshUseSpacetree
(
  const std::string& meshName ) const
{

  return _spacetreeNames.count(meshName) > 0;
}

const std::string& MeshConfiguration:: getSpacetreeName
(
  const std::string& meshName ) const
{
  assertion ( _spacetreeNames.count(meshName) > 0 );
  return _spacetreeNames.find(meshName)->second;
}

void MeshConfiguration:: addNeededMesh(
  const std::string& participant,
  const std::string& mesh)
{
  preciceTrace2 ( "addNeededMesh()", participant, mesh );
  if(_neededMeshes.count(participant)==0){
    std::vector<std::string> meshes;
    meshes.push_back(mesh);
    _neededMeshes.insert(std::pair<std::string,std::vector<std::string> >(participant, meshes));
  }
  else if(not utils::contained(mesh,_neededMeshes.find(participant)->second)){
    _neededMeshes.find(participant)->second.push_back(mesh);
  }
}

}} // namespace precice, mesh
