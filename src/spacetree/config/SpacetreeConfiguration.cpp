// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "SpacetreeConfiguration.hpp"
#include "spacetree/DynamicOctree.hpp"
#include "spacetree/StaticOctree.hpp"
//#include "spacetree/DynamicPeanotree2D.hpp"
//#include "spacetree/DynamicPeanotree3D.hpp"
#include "query/FindVoxelContent.hpp"
#include "utils/xml/XMLAttribute.hpp"
#include "utils/xml/ValidatorEquals.hpp"
#include "utils/xml/ValidatorOr.hpp"
#include "utils/Globals.hpp"
#include <string>

namespace precice {
namespace spacetree {

tarch::logging::Log SpacetreeConfiguration:: _log ( "precice::spacetree::SpacetreeConfiguration" );

//const std::string& SpacetreeConfiguration:: getTag()
//{
//  static std::string tag ( "spacetree" );
//  return tag;
//}

SpacetreeConfiguration:: SpacetreeConfiguration
(
  utils::XMLTag& parent )
:
  TAG("spacetree"),
  ATTR_NAME("name"),
  ATTR_TYPE("type"),
  VALUE_DYNAMIC_OCTREE("dynamic-octree"),
  VALUE_STATIC_OCTREE("static-octree"),
  VALUE_DYNAMIC_PEANOTREE2D("dynamic-peanotree-2D"),
  VALUE_DYNAMIC_PEANOTREE3D("dynamic-peanotree-3D"),
  _dimensions(0),
  _spacetrees()
  //_isValid(false)
{
  using namespace utils;
  XMLTag::Occurrence occ = XMLTag::OCCUR_ARBITRARY;
  std::list<XMLTag> tags;
  {
    XMLTag tag(*this, VALUE_DYNAMIC_OCTREE, occ, TAG);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_STATIC_OCTREE, occ, TAG);
    tags.push_back(tag);
  }
//  {
//    XMLTag tag(*this, VALUE_DYNAMIC_PEANOTREE2D, occ, TAG);
//    tags.push_back(tag);
//  }
//  {
//    XMLTag tag(*this, VALUE_DYNAMIC_PEANOTREE3D, occ, TAG);
//    tags.push_back(tag);
//  }

  utils::XMLAttribute<std::string> attrName ( ATTR_NAME );


//  utils::XMLAttribute<std::string> attrType ( ATTR_TYPE );
//  utils::ValidatorEquals<std::string> validDynamicOctree(VALUE_DYNAMIC_OCTREE);
//  utils::ValidatorEquals<std::string> validStaticOctree(VALUE_STATIC_OCTREE);
//  utils::ValidatorEquals<std::string> validDynamicPeanotree2D(VALUE_DYNAMIC_PEANOTREE2D);
//  utils::ValidatorEquals<std::string> validDynamicPeanotree3D(VALUE_DYNAMIC_PEANOTREE3D);
//  attrType.setValidator ( validDynamicOctree || validStaticOctree
//                          || validDynamicPeanotree2D || validDynamicPeanotree3D  );
//  tagSpacetree.addAttribute(attrType);


  utils::XMLAttribute<utils::DynVector> attrOffset("offset");
  attrOffset.setDefaultValue (utils::DynVector(3, 0.0));
  utils::XMLAttribute<utils::DynVector> attrHalflength("halflength");
  utils::XMLAttribute<double> attrMaxMeshwidth ("max-meshwidth");

  foreach (XMLTag& tag, tags){
    tag.addAttribute(attrName);
    tag.addAttribute(attrOffset);
    tag.addAttribute(attrHalflength);
    tag.addAttribute(attrMaxMeshwidth);
    parent.addSubtag(tag);
  }
}

void SpacetreeConfiguration:: setDimensions
(
  int dimensions )
{
  preciceTrace1("setDimensions()", dimensions);
  assertion1((dimensions == 2) || (dimensions == 3), dimensions);
  _dimensions = dimensions;
}

//bool SpacetreeConfiguration:: parseSubtag
//(
//  utils::XMLTag::XMLReader* xmlReader )
//{
//  utils::XMLTag tagSpacetree (TAG, utils::XMLTag::OCCUR_ONCE);
//
//  utils::XMLAttribute<std::string> attrName ( ATTR_NAME );
//  tagSpacetree.addAttribute ( attrName );
//
//  utils::XMLAttribute<std::string> attrType ( ATTR_TYPE );
//  utils::ValidatorEquals<std::string> validDynamicOctree(VALUE_DYNAMIC_OCTREE);
//  utils::ValidatorEquals<std::string> validStaticOctree(VALUE_STATIC_OCTREE);
//  utils::ValidatorEquals<std::string> validDynamicPeanotree2D(VALUE_DYNAMIC_PEANOTREE2D);
//  utils::ValidatorEquals<std::string> validDynamicPeanotree3D(VALUE_DYNAMIC_PEANOTREE3D);
//  attrType.setValidator ( validDynamicOctree || validStaticOctree
//                          || validDynamicPeanotree2D || validDynamicPeanotree3D  );
//  tagSpacetree.addAttribute(attrType);
//
//  if (_dimensions == 2){
//    using utils::Vector2D;
//    utils::XMLAttribute<Vector2D> attrOffset ("offset");
//    attrOffset.setDefaultValue (Vector2D(0.0));
//    tagSpacetree.addAttribute (attrOffset);
//
//    utils::XMLAttribute<Vector2D> attrHalflength ("halflength");
//    tagSpacetree.addAttribute (attrHalflength);
//  }
//  else {
//    using utils::Vector3D;
//    utils::XMLAttribute<Vector3D> attrOffset ("offset");
//    attrOffset.setDefaultValue (Vector3D(0.0));
//    tagSpacetree.addAttribute (attrOffset);
//
//    utils::XMLAttribute<Vector3D> attrHalflength ("halflength");
//    tagSpacetree.addAttribute (attrHalflength);
//  }
//
//  utils::XMLAttribute<double> attrMaxMeshwidth ("max-meshwidth");
//  tagSpacetree.addAttribute (attrMaxMeshwidth);
//
//  _isValid = tagSpacetree.parse ( xmlReader, *this );
//  return _isValid;
//}

const PtrSpacetree& SpacetreeConfiguration:: getSpacetree
(
  const std::string& name ) const
{
  //assertion ( _isValid );
  foreach ( const ConfiguredSpacetree& tree, _spacetrees ) {
    if ( tree.name == name ) {
      return tree.spacetree;
    }
  }
  preciceError ( "getSpacetree()", "A spacetree with name \"" << name
                 << "\" is not defined!" );
}

PtrSpacetree SpacetreeConfiguration:: getSpacetree
(
  const std::string&      type,
  const utils::DynVector& offset,
  const utils::DynVector& halflengths,
  double                  maxMeshwidth ) const
{
  preciceTrace4("getSpacetree()", type, offset, halflengths, maxMeshwidth);
  assertion(_dimensions != 0);
  Spacetree* spacetree = NULL;
  assertion2 ( offset.size() == halflengths.size(), offset.size(), halflengths.size() );
  bool equalHalflengths = true;
  equalHalflengths &= tarch::la::equals(halflengths(0), halflengths(1));
  if ( _dimensions == 3 ){
    equalHalflengths &= tarch::la::equals(halflengths(1), halflengths(2));
  }
  preciceCheck(equalHalflengths, "getSpacetree()", "All halflengths have to "
               << "be equal for a spacetree of type \"quad\"!");
  if (type == VALUE_DYNAMIC_OCTREE){
    spacetree = new DynamicOctree( offset, halflengths(0), maxMeshwidth);
  }
  else if (type == VALUE_STATIC_OCTREE){
    spacetree = new StaticOctree( offset, halflengths(0), maxMeshwidth);
  }
//  else if (type == VALUE_DYNAMIC_PEANOTREE2D){
//    spacetree = new DynamicPeanotree2D( offset, halflengths(0), maxMeshwidth);
//  }
//  else if (type == VALUE_DYNAMIC_PEANOTREE3D){
//    spacetree = new DynamicPeanotree3D( offset, halflengths(0), maxMeshwidth);
//  }
  else {
    preciceError("getSpacetree()", "Unknown spacetree type \"" << type << "\"!");
  }
  assertion(spacetree != NULL);
  return PtrSpacetree(spacetree);
}

void SpacetreeConfiguration:: xmlTagCallback
(
  utils::XMLTag& tag )
{
  preciceTrace1 ( "xmlTagCallback()", tag.getName() );
  if (tag.getNamespace() == TAG){
    assertion(_dimensions != 0);
    std::string name = tag.getStringAttributeValue(ATTR_NAME);
    std::string type = tag.getName();
    utils::DynVector offset(_dimensions);
    utils::DynVector halflengths(_dimensions);
    offset = tag.getDynVectorAttributeValue("offset", _dimensions);
    halflengths = tag.getDynVectorAttributeValue("halflength", _dimensions);
    double maxMeshwidth = tag.getDoubleAttributeValue("max-meshwidth");
    foreach (const ConfiguredSpacetree& tree, _spacetrees){
      if (tree.name == name){
        std::ostringstream stream;
        stream << "Spacetree with " << "name \"" << name << "\" is defined twice";
        throw stream.str();
      }
    }
    ConfiguredSpacetree tree;
    tree.name = name;
    tree.spacetree = getSpacetree(type, offset, halflengths, maxMeshwidth);
    _spacetrees.push_back(tree);
  }
}

}} // namespace precice, spacetree
