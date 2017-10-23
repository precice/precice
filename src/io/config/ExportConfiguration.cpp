#include "ExportConfiguration.hpp"
#include "io/Export.hpp"
#include "utils/Globals.hpp"
#include "utils/Helpers.hpp"
#include "utils/xml/XMLTag.hpp"
#include "utils/xml/XMLAttribute.hpp"
#include "utils/xml/ValidatorEquals.hpp"
#include "utils/xml/ValidatorOr.hpp"

namespace precice {
namespace io {

logging::Logger ExportConfiguration:: _log("io::ExportConfiguration");

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
    tag.addAttribute(attrEveryIteration);
    parent.addSubtag(tag);
  }
}

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
    context.plotNormals = tag.getBooleanAttributeValue(ATTR_NORMALS);
    context.everyIteration = tag.getBooleanAttributeValue(ATTR_EVERY_ITERATION);
    context.type = tag.getName();
    if ((context.timestepInterval == -1) &&  context.triggerSolverPlot){
      std::string error = "Attribute timestep interval has to be set when ";
      error += "trigger-solver is activated";
      throw error;
    }
    _contexts.push_back(context);
  }
}


}} // namespace precice, io

