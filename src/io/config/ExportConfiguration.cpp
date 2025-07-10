#include "ExportConfiguration.hpp"
#include "xml/ConfigParser.hpp"
#include "xml/XMLAttribute.hpp"
#include "xml/XMLTag.hpp"

namespace precice::io {

ExportConfiguration::ExportConfiguration(xml::XMLTag &parent)
{
  using namespace xml;
  std::string        doc;
  std::list<XMLTag>  tags;
  XMLTag::Occurrence occ = XMLTag::OCCUR_ARBITRARY;
  {
    XMLTag tag(*this, VALUE_VTK, occ, TAG);
    tag.setDocumentation("Exports meshes to VTK legacy format files. Parallel participants will use the VTU exporter instead.");
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_VTU, occ, TAG);
    tag.setDocumentation("Exports meshes to VTU files in serial or PVTU files with VTU piece files in parallel.");
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_VTP, occ, TAG);
    tag.setDocumentation("Exports meshes to VTP files in serial or PVTP files with VTP piece files in parallel.");
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_CSV, occ, TAG);
    tag.setDocumentation("Exports vertex coordinates and data to CSV files.");
    tags.push_back(tag);
  }

  auto attrLocation = XMLAttribute<std::string>(ATTR_LOCATION, ".")
                          .setDocumentation("Directory to export the files to.");

  auto attrEveryNTimeWindows = makeXMLAttribute(ATTR_EVERY_N_TIME_WINDOWS, 1)
                                   .setDocumentation("preCICE does an export every X time windows. Choose -1 for no exports.");

  auto attrEveryIteration = makeXMLAttribute(ATTR_EVERY_ITERATION, false)
                                .setDocumentation("Exports in every coupling (sub)iteration. For debug purposes.");

  auto attrUpdateSeries = makeXMLAttribute(ATTR_UPDATE_SERIES, false)
                              .setDocumentation("Update the series file after every export instead of at the end of the simulation.");

  for (XMLTag &tag : tags) {
    tag.addAttribute(attrLocation);
    tag.addAttribute(attrEveryNTimeWindows);
    tag.addAttribute(attrEveryIteration);
    tag.addAttribute(attrUpdateSeries);
    parent.addSubtag(tag);
  }
}

void ExportConfiguration::xmlTagCallback(
    const xml::ConfigurationContext &context,
    xml::XMLTag                     &tag)
{
  if (tag.getNamespace() == TAG) {
    ExportContext econtext;
    econtext.location          = tag.getStringAttributeValue(ATTR_LOCATION);
    econtext.everyNTimeWindows = tag.getIntAttributeValue(ATTR_EVERY_N_TIME_WINDOWS);
    econtext.everyIteration    = tag.getBooleanAttributeValue(ATTR_EVERY_ITERATION);
    econtext.updateSeries      = tag.getBooleanAttributeValue(ATTR_UPDATE_SERIES);
    econtext.type              = tag.getName();
    _contexts.push_back(econtext);
  }
}

} // namespace precice::io
