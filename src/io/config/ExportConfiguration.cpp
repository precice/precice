#include "ExportConfiguration.hpp"
#include "io/Export.hpp"
#include "xml/XMLAttribute.hpp"
#include "xml/XMLTag.hpp"

namespace precice {
namespace io {

ExportConfiguration::ExportConfiguration(xml::XMLTag &parent)
{
  using namespace xml;
  std::string        doc;
  std::list<XMLTag>  tags;
  XMLTag::Occurrence occ = XMLTag::OCCUR_ARBITRARY;
  {
    XMLTag tag(*this, VALUE_VTK, occ, TAG);
    tag.setDocumentation("Exports meshes to VTK text files.");
    tags.push_back(tag);
  }

  auto attrLocation = XMLAttribute<std::string>(ATTR_LOCATION, "")
                          .setDocumentation("Directory to export the files to.");

  auto attrTimestepInterval = makeXMLAttribute(ATTR_TIMESTEP_INTERVAL, 1)
                                  .setDocumentation("preCICE timestep interval for export of files. Choose -1 for no exports.");

  auto attrTriggerSolver = makeXMLAttribute(ATTR_TRIGGER_SOLVER, false)
                               .setDocumentation(
                                   std::string("If set to on/yes, an action requirement is set for the participant ") +
                                   "with frequency defined by attribute " + ATTR_TIMESTEP_INTERVAL + ".");

  auto attrNormals = makeXMLAttribute(ATTR_NORMALS, true)
                         .setDocumentation("If set to on/yes, mesh normals (if available) are added to the export.");

  auto attrEveryIteration = makeXMLAttribute(ATTR_EVERY_ITERATION, false)
                                .setDocumentation("Exports in every coupling (sub)iteration. For debug purposes.");

  for (XMLTag &tag : tags) {
    tag.addAttribute(attrLocation);
    tag.addAttribute(attrTimestepInterval);
    tag.addAttribute(attrTriggerSolver);
    tag.addAttribute(attrNormals);
    tag.addAttribute(attrEveryIteration);
    parent.addSubtag(tag);
  }
}

void ExportConfiguration::xmlTagCallback(
    const xml::ConfigurationContext& context,
    xml::XMLTag &tag)
{
  if (tag.getNamespace() == TAG) {
    ExportContext context;
    context.location          = tag.getStringAttributeValue(ATTR_LOCATION);
    context.triggerSolverPlot = tag.getBooleanAttributeValue(ATTR_TRIGGER_SOLVER);
    context.timestepInterval  = tag.getIntAttributeValue(ATTR_TIMESTEP_INTERVAL);
    context.plotNormals       = tag.getBooleanAttributeValue(ATTR_NORMALS);
    context.everyIteration    = tag.getBooleanAttributeValue(ATTR_EVERY_ITERATION);
    context.type              = tag.getName();
    _contexts.push_back(context);
  }
}

} // namespace io
} // namespace precice
