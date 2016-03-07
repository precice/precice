// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "ExportConfiguration.hpp"
#include "io/ExportVTK.hpp"
#include "io/ExportVRML.hpp"
#include "io/Export.hpp"
#include "utils/Globals.hpp"
#include "utils/Helpers.hpp"
#include "utils/xml/XMLTag.hpp"
#include "utils/xml/XMLAttribute.hpp"
#include "utils/xml/ValidatorEquals.hpp"
#include "utils/xml/ValidatorOr.hpp"

namespace precice {
namespace io {

tarch::logging::Log ExportConfiguration:: _log("precice::io::ExportConfiguration");

//const std::string& ExportConfiguration:: getTag()
//{
//  static std::string tag("export");
//  return tag;
//}

ExportConfiguration:: ExportConfiguration
(
  utils::XMLTag& parent )
:
  TAG("export"),
  ATTR_LOCATION ( "directory" ),
  ATTR_TYPE ( "type" ),
  ATTR_AUTO ( "auto" ),
  VALUE_VTK ( "vtk" ),
  VALUE_VRML ( "vrml" ),
  ATTR_TIMESTEP_INTERVAL ( "timestep-interval" ),
  ATTR_NEIGHBORS ( "neighbors" ),
  ATTR_TRIGGER_SOLVER ( "trigger-solver" ),
  ATTR_NORMALS ( "normals" ),
  ATTR_SPACETREE ( "spacetree" ),
  ATTR_EVERY_ITERATION("every-iteration"),
  //_isValid ( false ),
  _contexts()
{
  using namespace utils;
  std::string doc;
  std::list<XMLTag> tags;
  XMLTag::Occurrence occ = XMLTag::OCCUR_ARBITRARY;
  {
    XMLTag tag(*this, VALUE_VTK, occ, TAG);
    tag.setDocumentation("Exports meshes to VTK text files.");
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_VRML, occ, TAG);
    tag.setDocumentation("Exports meshes to VRML 1.0 text files.");
    tags.push_back(tag);
  }

  XMLAttribute<std::string> attrLocation(ATTR_LOCATION);
  attrLocation.setDocumentation("Directory to export the files to.");
  attrLocation.setDefaultValue("");

  XMLAttribute<int> attrTimestepInterval(ATTR_TIMESTEP_INTERVAL);
  doc = "preCICE timestep interval for export of files. If set to -1 (default), ";
  doc += "files are only exported when a solver requests the export.";
  attrTimestepInterval.setDocumentation(doc);
  attrTimestepInterval.setDefaultValue(-1);

//  XMLAttribute<bool> attrNeighbors(ATTR_NEIGHBORS);
//  attrNeighbors.setDefaultValue(false);

  XMLAttribute<bool> attrTriggerSolver(ATTR_TRIGGER_SOLVER);
  doc = "If set to on/yes, an action requirement is set for the participant ";
  doc += "with frequency defined by attribute " + ATTR_TIMESTEP_INTERVAL + ".";
  attrTriggerSolver.setDocumentation(doc);
  attrTriggerSolver.setDefaultValue(false);

  XMLAttribute<bool> attrNormals(ATTR_NORMALS);
  doc = "If set to on/yes, mesh normals are added to the export.";
  attrNormals.setDocumentation(doc);
  attrNormals.setDefaultValue(true);

  XMLAttribute<bool> attrSpacetree(ATTR_SPACETREE);
  doc = "If set to on/yes, spacetrees used by a mesh are also exported.";
  attrSpacetree.setDocumentation(doc);
  attrSpacetree.setDefaultValue(true);

  XMLAttribute<bool> attrEveryIteration(ATTR_EVERY_ITERATION);
  doc = "Exports in every coupling (sub)iteration. For debug purposes.";
  attrEveryIteration.setDocumentation(doc);
  attrEveryIteration.setDefaultValue(false);

  for (XMLTag& tag : tags){
    tag.addAttribute(attrLocation);
    tag.addAttribute(attrTimestepInterval);
    //tag.addAttribute(attrNeighbors);
    tag.addAttribute(attrTriggerSolver);
    tag.addAttribute(attrNormals);
    tag.addAttribute(attrSpacetree);
    tag.addAttribute(attrEveryIteration);
    parent.addSubtag(tag);
  }
}

//bool ExportConfiguration:: parseSubtag
//(
//  utils::XMLTag::XMLReader* xmlReader )
//{
//  using utils::XMLTag;
//  using utils::XMLAttribute;
//  using utils::ValidatorEquals;
//  XMLTag tag ( TAG, XMLTag::OCCUR_ARBITRARY );
//
//  XMLAttribute<std::string> attrLocation ( ATTR_LOCATION );
//  attrLocation.setDefaultValue ("");
//  tag.addAttribute ( attrLocation );
//
//  XMLAttribute<std::string> attrType ( ATTR_TYPE );
//  ValidatorEquals<std::string> validTypeVTK ( VALUE_VTK );
//  ValidatorEquals<std::string> validTypeSocket ( VALUE_VRML );
//  attrType.setValidator ( validTypeVTK || validTypeSocket );
//  attrType.setDefaultValue ( VALUE_VTK );
//  tag.addAttribute ( attrType );
//
//  XMLAttribute<int> attrTimestepInterval ( ATTR_TIMESTEP_INTERVAL );
//  attrTimestepInterval.setDefaultValue ( -1 );
//  tag.addAttribute ( attrTimestepInterval );
//
//  XMLAttribute<bool> attrNeighbors ( ATTR_NEIGHBORS );
//  attrNeighbors.setDefaultValue ( false );
//  tag.addAttribute ( attrNeighbors );
//
//  XMLAttribute<bool> attrTriggerSolver ( ATTR_TRIGGER_SOLVER );
//  attrTriggerSolver.setDefaultValue ( false );
//  tag.addAttribute ( attrTriggerSolver );
//
//  XMLAttribute<bool> attrNormals ( ATTR_NORMALS );
//  attrNormals.setDefaultValue ( false );
//  tag.addAttribute ( attrNormals );
//
//  XMLAttribute<bool> attrSpacetree ( ATTR_SPACETREE );
//  attrSpacetree.setDefaultValue(true);
//  tag.addAttribute(attrSpacetree);
//
//  //utils::XMLAttribute<bool> attrAuto ( ATTR_AUTO );
//  //attrAuto.setDefaultValue ( true );
//  //tagExport.addAttribute ( attrAuto );
//
//  _context = ExportContext();
//  _isValid = tag.parse ( xmlReader, *this );
//
//  return _isValid;
//}

void ExportConfiguration:: xmlTagCallback
(
  utils::XMLTag& tag )
{
  if ( tag.getNamespace() == TAG ){
    ExportContext context;
    context.location = tag.getStringAttributeValue(ATTR_LOCATION);
    //context.plotNeighbors = tag.getBooleanAttributeValue(ATTR_NEIGHBORS);
    context.triggerSolverPlot =  tag.getBooleanAttributeValue(ATTR_TRIGGER_SOLVER);
    context.timestepInterval = tag.getIntAttributeValue(ATTR_TIMESTEP_INTERVAL);
    bool plotNormals = tag.getBooleanAttributeValue(ATTR_NORMALS);
    context.exportSpacetree = tag.getBooleanAttributeValue(ATTR_SPACETREE);
    context.everyIteration = tag.getBooleanAttributeValue(ATTR_EVERY_ITERATION);
    std::string type = tag.getName();
    if ((context.timestepInterval == -1) &&  context.triggerSolverPlot){
      std::string error = "Attribute timestep interval has to be set when ";
      error += "trigger-solver is activated";
      throw error;
    }
    PtrExport exporter;
    if (type == VALUE_VTK){
      exporter = PtrExport(new ExportVTK(plotNormals));
    }
    else if (type == VALUE_VRML){
      exporter = PtrExport (new ExportVRML(plotNormals));
    }
    else {
      preciceError("xmlTagCallback()", "Unknown export type!");
    }
    context.exporter = exporter;
    _contexts.push_back(context);
  }
}

}} // namespace precice, io

