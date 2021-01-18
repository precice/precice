#include "ExportConfiguration.hpp"
#include "xml/ConfigParser.hpp"
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

  auto attrEveryNTimeWindows = makeXMLAttribute(ATTR_EVERY_N_TIME_WINDOWS, 1)
                                   .setDocumentation("preCICE does an export every X time windows. Choose -1 for no exports.");

  auto attrNormals = makeXMLAttribute(ATTR_NORMALS, true)
                         .setDocumentation("If set to on/yes, mesh normals (if available) are added to the export.");

  auto attrEveryIteration = makeXMLAttribute(ATTR_EVERY_ITERATION, false)
                                .setDocumentation("Exports in every coupling (sub)iteration. For debug purposes.");

  for (XMLTag &tag : tags) {
    tag.addAttribute(attrLocation);
    tag.addAttribute(attrEveryNTimeWindows);
    tag.addAttribute(attrNormals);
    tag.addAttribute(attrEveryIteration);
    parent.addSubtag(tag);
  }
}

void ExportConfiguration::xmlTagCallback(
    const xml::ConfigurationContext &context,
    xml::XMLTag &                    tag)
{
  if (tag.getNamespace() == TAG) {
    ExportContext econtext;
    econtext.location          = tag.getStringAttributeValue(ATTR_LOCATION);
    econtext.everyNTimeWindows = tag.getIntAttributeValue(ATTR_EVERY_N_TIME_WINDOWS);
    econtext.plotNormals       = tag.getBooleanAttributeValue(ATTR_NORMALS);
    econtext.everyIteration    = tag.getBooleanAttributeValue(ATTR_EVERY_ITERATION);
    econtext.type              = tag.getName();
    _contexts.push_back(econtext);
  }
}

} // namespace io
} // namespace precice
