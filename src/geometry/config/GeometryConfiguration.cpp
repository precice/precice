#include "GeometryConfiguration.hpp"
#include "geometry/Sphere.hpp"
#include "geometry/Bubble.hpp"
#include "geometry/Cuboid.hpp"
#include "geometry/DriftRatchet.hpp"
#include "geometry/ImportGeometry.hpp"
#include "geometry/SharedPointer.hpp"
#include "mesh/config/DataConfiguration.hpp"
#include "mesh/config/MeshConfiguration.hpp"
#include "mesh/Mesh.hpp"
#include "com/config/CommunicationConfiguration.hpp"
#include "utils/xml/XMLAttribute.hpp"
#include "utils/xml/ValidatorEquals.hpp"
#include "utils/xml/ValidatorOr.hpp"
#include "utils/Globals.hpp"
#include <iostream>

namespace precice {
namespace geometry {

logging::Logger GeometryConfiguration:: _log ("precice::geometry::GeometryConfiguration");

GeometryConfiguration:: GeometryConfiguration
(
  utils::XMLTag&             parent,
  mesh::PtrMeshConfiguration meshConfig )
:
  TAG("geometry"),
  TAG_LENGTH("length"),
  TAG_DISCRETIZATION_WIDTH("discretization-width"),
  TAG_OFFSET("offset"),
  TAG_PORES("pores"),
  TAG_RADIUS("radius"),
  TAG_DEFORMATION("deformation"),
  TAG_FILENAME("filename"),
  TAG_FILETYPE("filetype"),
  TAG_PROVIDER("provider"),
  TAG_RECEIVER("receiver"),
  ATTR_MESH("of-mesh"),
  ATTR_TYPE("type"),
  ATTR_VALUE("value"),
  ATTR_NAME("name"),
  VALUE_BUILTIN_CUBOID("builtin-cuboid"),
  VALUE_BUILTIN_SPHERE("builtin-sphere"),
  VALUE_BUILTIN_BUBBLE("builtin-bubble"),
  VALUE_BUILTIN_DRATCHET("builtin-dratchet"),
  VALUE_BUILTIN_FACE("builtin-face"),
  VALUE_IMPORT("import"),
  VALUE_NONE("none"),
  VALUE_AUTO("auto"),
  _meshConfig(meshConfig),
  _dimensions(0),
  _readData(_dimensions),
  _geometries()
  //_isValid(false)
{
  using namespace utils;
  std::string doc;
  std::list<XMLTag> tags;
  XMLTag::Occurrence occ = XMLTag::OCCUR_ARBITRARY;
  {
    XMLTag tag(*this, VALUE_BUILTIN_CUBOID, occ, TAG);
    doc = "Rectangle in 2d, hexahedron in 3d.";
    tag.setDocumentation(doc);
    addTypeSpecificAttributes(VALUE_BUILTIN_CUBOID, tag);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_BUILTIN_SPHERE, occ, TAG);
    doc = "Circle in 2d, sphere in 3d.";
    tag.setDocumentation(doc);
    addTypeSpecificAttributes(VALUE_BUILTIN_SPHERE, tag);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_BUILTIN_BUBBLE, occ, TAG);
    doc = "Bubble-like geometry.";
    tag.setDocumentation(doc);
    addTypeSpecificAttributes(VALUE_BUILTIN_BUBBLE, tag);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_BUILTIN_DRATCHET, occ, TAG);
    doc = "In 2d a channel with periodically varying height, in 3d a cylinder ";
    doc += "with periodically varying diameter.";
    tag.setDocumentation(doc);
    addTypeSpecificAttributes(VALUE_BUILTIN_DRATCHET, tag);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_IMPORT, occ, TAG);
    doc = "Imports a geometry from a VRML 1.0 file.";
    tag.setDocumentation(doc);
    addTypeSpecificAttributes(VALUE_IMPORT, tag);
    tags.push_back(tag);
  }

  XMLAttribute<std::string> attrMesh(ATTR_MESH);
  attrMesh.setDocumentation("Name of mesh to build the geometry from.");

//  XMLAttribute<std::string> attrGeometryType(ATTR_TYPE);
//  ValidatorEquals<std::string> validCuboidType(VALUE_BUILTIN_CUBOID);
//  ValidatorEquals<std::string> validSphereType(VALUE_BUILTIN_SPHERE);
//  ValidatorEquals<std::string> validBubbleType(VALUE_BUILTIN_BUBBLE);
//  ValidatorEquals<std::string> validDratchetType(VALUE_BUILTIN_DRATCHET);
//  ValidatorEquals<std::string> validImportGeometry(VALUE_IMPORT);
//  attrGeometryType.setValidator(validCuboidType || validSphereType
//      || validBubbleType || validDratchetType || validImportGeometry);
//  _xmlTag.addAttribute(attrGeometryType);

  XMLTag tagOffset(*this, TAG_OFFSET, XMLTag::OCCUR_NOT_OR_ONCE);
  XMLAttribute<utils::DynVector> attrValue(ATTR_VALUE);
  tagOffset.addAttribute(attrValue);

  for (XMLTag& tag : tags) {
    tag.addAttribute(attrMesh);
    tag.addSubtag(tagOffset);
    parent.addSubtag(tag);
  }
}

void GeometryConfiguration:: setDimensions
(
  int dimensions )
{
  preciceTrace("setDimensions()", dimensions);
  assertion((dimensions == 2) || (dimensions == 3), dimensions);
  _dimensions = dimensions;
}

//bool GeometryConfiguration:: parseSubtag
//(
//  utils::XMLTag::XMLReader* xmlReader )
//{
//  preciceTrace ( "parseSubtag()" );
////  using utils::XMLTag;
////  using utils::XMLAttribute;
////  using utils::ValidatorEquals;
//
//  //if ( _xmlTag.isConfigured() ) {  // This is not the first geometry tag
//  //  _xmlTag.clear ();
//  //}
//
////  XMLAttribute<std::string> attrMesh ( ATTR_MESH );
////  _xmlTag.addAttribute ( attrMesh );
////
////  XMLAttribute<std::string> attrGeometryType ( ATTR_TYPE );
////  ValidatorEquals<std::string> validCuboidType ( VALUE_BUILTIN_CUBOID );
////  ValidatorEquals<std::string> validSphereType ( VALUE_BUILTIN_SPHERE );
////  ValidatorEquals<std::string> validBubbleType ( VALUE_BUILTIN_BUBBLE );
////  ValidatorEquals<std::string> validDratchetType ( VALUE_BUILTIN_DRATCHET );
////  ValidatorEquals<std::string> validImportGeometry ( VALUE_IMPORT );
////  attrGeometryType.setValidator ( validCuboidType || validSphereType
////      || validBubbleType || validDratchetType || validImportGeometry );
////  _xmlTag.addAttribute ( attrGeometryType );
////
////  XMLTag tagOffset ( TAG_OFFSET, XMLTag::OCCUR_NOT_OR_ONCE );
////  if (_dimensions == 2){
////    XMLAttribute<utils::Vector2D> attrValue ( ATTR_VALUE );
////    tagOffset.addAttribute ( attrValue );
////  }
////  else {
////    XMLAttribute<utils::Vector3D> attrValue ( ATTR_VALUE );
////    tagOffset.addAttribute ( attrValue );
////  }
////  _xmlTag.addSubtag ( tagOffset );
//
//  _isValid = _tag.parse(xmlReader);
//  return _isValid;
//}

void GeometryConfiguration:: xmlTagCallback
(
  utils::XMLTag& tag )
{
  preciceTrace("xmlTagCallback()", tag.getFullName());
  if (tag.getNamespace() == TAG){
    assertion(_dimensions != 0);
    _readData = ReadData(_dimensions);
    _readData.noReadData = false;
    _readData.type = tag.getName();
    _readData.mesh = tag.getStringAttributeValue(ATTR_MESH);
    //addTypeSpecificAttributes ( _readData.type, tag );
  }
  else if (tag.getName() == TAG_OFFSET){
    assertion(not _readData.noReadData);
    assertion(_dimensions != 0);
    _readData.offset = tag.getDynVectorAttributeValue(ATTR_VALUE, _dimensions);
  }
  else if (tag.getName() == TAG_RADIUS){
    assertion(not _readData.noReadData);
    _readData.radius = tag.getDoubleAttributeValue(ATTR_VALUE);
  }
  else if (tag.getName() == TAG_DEFORMATION){
    assertion(not _readData.noReadData);
    _readData.deformation = tag.getDoubleAttributeValue(ATTR_VALUE);
  }
  else if (tag.getName() == TAG_DISCRETIZATION_WIDTH){
    assertion(not _readData.noReadData);
    _readData.discretizationWidth = tag.getDoubleAttributeValue(ATTR_VALUE);
  }
  else if (tag.getName() == TAG_LENGTH){
    assertion(not _readData.noReadData);
    if (_readData.type == VALUE_BUILTIN_DRATCHET){
      _readData.scalarLength = tag.getDoubleAttributeValue(ATTR_VALUE);
    }
    else if (_readData.type == VALUE_BUILTIN_CUBOID){
      assertion(_dimensions != 0);
      _readData.length = tag.getDynVectorAttributeValue(ATTR_VALUE, _dimensions);
    }
    else {
      assertion(false);
    }
  }
  else if (tag.getName() == TAG_PORES){
    assertion(not _readData.noReadData);
    assertion(_readData.type == VALUE_BUILTIN_DRATCHET);
    _readData.pores = tag.getDoubleAttributeValue(ATTR_VALUE);
  }
  else if (tag.getName() == TAG_FILENAME){
    assertion(not _readData.noReadData);
    _readData.filename = tag.getStringAttributeValue(ATTR_VALUE);
  }
  else if (tag.getName() == TAG_FILETYPE){
    assertion(not _readData.noReadData);
    _readData.filetype = tag.getStringAttributeValue(ATTR_VALUE);
  }
}

void GeometryConfiguration:: xmlEndTagCallback
(
  utils::XMLTag& tag )
{
  preciceTrace ( "xmlEndTagCallback()", tag.getFullName() );
  if (tag.getNamespace() == TAG ) {
    assertion(not _readData.noReadData);
    assertion(_readData.type == tag.getName(), _readData.type, tag.getName());
    if (tag.getName() == VALUE_BUILTIN_CUBOID){
      addCuboid();
    }
    else if (tag.getName() == VALUE_BUILTIN_DRATCHET){
      addDriftRatchet();
    }
    else if (tag.getName() == VALUE_BUILTIN_SPHERE){
      addSphere();
    }
    else if (tag.getName() == VALUE_BUILTIN_BUBBLE){
      addBubble();
    }
    else if (tag.getName() == VALUE_IMPORT){
      addImportGeometry();
    }
    else {
      assertion(false);
    }
  }
}

void GeometryConfiguration:: addTypeSpecificAttributes
(
  const std::string& type,
  utils::XMLTag&     tag )
{
  preciceTrace("addTypeSpecificAttributes()", type);
  using utils::XMLTag;
  using utils::XMLAttribute;
  using utils::ValidatorEquals;

  XMLTag tagDiscretizationWidth(*this, TAG_DISCRETIZATION_WIDTH, XMLTag::OCCUR_ONCE);
  XMLAttribute<double> attrDoubleValue(ATTR_VALUE);
  tagDiscretizationWidth.addAttribute(attrDoubleValue);

  if (type == VALUE_BUILTIN_CUBOID){
    tag.addSubtag(tagDiscretizationWidth);

    XMLTag tagLength(*this, TAG_LENGTH, XMLTag::OCCUR_ONCE);
    XMLAttribute<utils::DynVector> attrValue(ATTR_VALUE);
    tagLength.addAttribute(attrValue);
    tag.addSubtag(tagLength);
  }
  else if (type == VALUE_BUILTIN_DRATCHET){
    tag.addSubtag(tagDiscretizationWidth);

    XMLTag tagLength (*this, TAG_LENGTH, XMLTag::OCCUR_ONCE);
    XMLAttribute<double> attrValue(ATTR_VALUE);
    //ValidatorGreaterThan<double> validValue (0.0);
    //attrValue.setValidator ( validValue );
    tagLength.addAttribute(attrValue);
    tag.addSubtag(tagLength);

    XMLTag tagPores(*this, TAG_PORES, XMLTag::OCCUR_ONCE);
    tagPores.addAttribute(attrValue);
    tag.addSubtag(tagPores);

    XMLTag tagRadius(*this, TAG_RADIUS, XMLTag::OCCUR_ONCE);
    tagRadius.addAttribute(attrValue);
    tag.addSubtag(tagRadius);
  }
  else if (type == VALUE_BUILTIN_SPHERE){
    tag.addSubtag(tagDiscretizationWidth);

    XMLTag tagRadius(*this, TAG_RADIUS, XMLTag::OCCUR_ONCE);
    XMLAttribute<double> attrValue(ATTR_VALUE);
    //ValidatorGreaterThan<double> validValue ( 0.0 );
    //attrValue.setValidator ( validValue );
    tagRadius.addAttribute(attrValue);
    tag.addSubtag(tagRadius);
  }
  else if (type == VALUE_BUILTIN_BUBBLE){
    tag.addSubtag(tagDiscretizationWidth);

    XMLTag tagRadius(*this, TAG_RADIUS, XMLTag::OCCUR_ONCE);
    XMLAttribute<double> attrValue(ATTR_VALUE);
    //ValidatorGreaterThan<double> validValue ( 0.0 );
    //attrValue.setValidator ( validValue );
    tagRadius.addAttribute(attrValue);
    tag.addSubtag(tagRadius);

    XMLTag tagDefo(*this, TAG_DEFORMATION, XMLTag::OCCUR_NOT_OR_ONCE);
    XMLAttribute<double> attrValue2(ATTR_VALUE);
    //ValidatorGreaterThan<double> validValue2 ( 0.0 );
    //attrValue2.setValidator ( validValue2 );
    tagDefo.addAttribute(attrValue2);
    tag.addSubtag(tagDefo);
  }
  else if (type == VALUE_IMPORT){
    XMLTag tagFilename(*this, TAG_FILENAME, XMLTag::OCCUR_ONCE);
    XMLAttribute<std::string> attrValue(ATTR_VALUE);
    tagFilename.addAttribute(attrValue);
    tag.addSubtag(tagFilename);

    XMLTag tagFiletype(*this, TAG_FILETYPE, XMLTag::OCCUR_ONCE);
    ValidatorEquals<std::string> validFiletype("vrml");
    attrValue.setValidator(validFiletype);
    tagFiletype.addAttribute(attrValue);
    tag.addSubtag(tagFiletype);
  }
  else {
    assertion (false);
  }
}

bool GeometryConfiguration:: addCuboid()
{
  assertion(_dimensions != 0);
  assertion ( _readData.mesh != std::string("") );
  assertion ( _readData.discretizationWidth > 0.0 );
  assertion ( not math::equals(_readData.length, Eigen::VectorXd::Zero(_dimensions)) );
  checkMeshName ( _readData.mesh );
  Cuboid* cuboidGeometry = new Cuboid ( _readData.offset,
      _readData.discretizationWidth,  _readData.length );
  _geometries.push_back ( PtrGeometry(cuboidGeometry) );
  _meshNames.push_back ( _readData.mesh );
//  _readData = ReadData (_dimensions);
  return true;
}

bool GeometryConfiguration:: addDriftRatchet()
{
  assertion(_dimensions != 0);
  assertion(_readData.mesh != std::string(""));
  assertion(_readData.scalarLength > 0.0);
  assertion(_readData.pores > 0.0);
  assertion(_readData.discretizationWidth > 0.0);
  assertion(_readData.radius > 0.0);

  checkMeshName(_readData.mesh);

  DriftRatchet* ratchetGeometry = new DriftRatchet (
      _readData.offset,
      _readData.discretizationWidth,
      _readData.radius,
      DriftRatchet::getDefaultMinRadius ( _readData.radius ),
      DriftRatchet::getDefaultShapeParameter (),
      _readData.scalarLength,
      _readData.pores,
      0, 0, 0 );

  _geometries.push_back(PtrGeometry(ratchetGeometry));
  _meshNames.push_back(_readData.mesh);
//  _readData = ReadData(_dimensions);
  return true;
}

void GeometryConfiguration:: addSphere()
{
  assertion(_dimensions != 0);
  assertion(_readData.mesh != std::string(""));
  assertion(_readData.radius > 0.0);
  assertion(_readData.discretizationWidth > 0.0);

  checkMeshName(_readData.mesh);

  Sphere* sphereGeometry = new Sphere (
      _readData.offset, _readData.discretizationWidth, _readData.radius );
  _geometries.push_back ( PtrGeometry(sphereGeometry) );
  _meshNames.push_back ( _readData.mesh );
//  _readData = ReadData (_dimensions);
}

void GeometryConfiguration:: addBubble()
{
  assertion(_dimensions != 0);
  assertion(_readData.mesh != std::string(""));
  assertion(_readData.radius > 0.0);
  assertion(_readData.deformation > -4/3.0);
  assertion(_readData.deformation <  4/3.0);
  assertion(_readData.discretizationWidth > 0.0);

  checkMeshName(_readData.mesh);

  Bubble* bubbleGeometry = new Bubble (
      _readData.offset, _readData.discretizationWidth,
      _readData.radius, _readData.deformation );
  _geometries.push_back(PtrGeometry(bubbleGeometry));
  _meshNames.push_back(_readData.mesh);
//  _readData = ReadData(_dimensions);
}

void GeometryConfiguration:: addImportGeometry()
{
  assertion(_dimensions != 0);
  assertion(_readData.mesh != std::string(""));
  assertion(_readData.filename != std::string(""));
  assertion(_readData.filetype != std::string(""));

  checkMeshName(_readData.mesh);

  ImportGeometry* importGeometry = nullptr;
  if (_readData.filetype == std::string("vrml")){
    bool importCheckpoint = false; // Only mesh topology is imported.
    importGeometry = new ImportGeometry(_readData.offset, _readData.filename,
                                        ImportGeometry::VRML_1_FILE, importCheckpoint, true);
  }
  assertion(importGeometry != nullptr);
  _geometries.push_back(PtrGeometry(importGeometry));
  _meshNames.push_back(_readData.mesh);
//  _readData = ReadData(_dimensions);
}

void GeometryConfiguration:: checkMeshName
(
   const std::string& meshName )
{
  bool found = false;
  for (const mesh::PtrMesh mesh : _meshConfig->meshes()) {
    if ( mesh->getName() == meshName ) {
      found = true;
      break;
    }
  }
  preciceCheck ( found, "checkMesh()", "Geometry references mesh \""
      << meshName << "\", which is not given in configuration!" );
}


//bool GeometryConfiguration:: isValid() const
//{
//  return _isValid;
//}

PtrGeometry GeometryConfiguration:: getGeometry
(
  const std::string& meshName )
{
  for ( size_t i=0; i < _meshNames.size(); i++ ) {
    if ( _meshNames[i] == meshName ) {
      return _geometries[i];
    }
  }
  return PtrGeometry();
}

void GeometryConfiguration:: addGeometry
(
  PtrGeometry        geometry,
  const std::string& meshName)
{
  _geometries.push_back ( geometry );
  preciceCheck ( ! utils::contained(meshName, _meshNames),
      "addGeometry()", "Mesh name is not unique!" );
  _meshNames.push_back ( meshName );
}

}} // namespace precice, geometry
