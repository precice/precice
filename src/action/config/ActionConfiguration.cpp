#include "ActionConfiguration.hpp"
#include <algorithm>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <utility>

#include "action/PythonAction.hpp"
#include "action/RecorderAction.hpp"
#include "action/ScaleByAreaAction.hpp"
#include "action/SummationAction.hpp"
#include "logging/LogMacros.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/config/MeshConfiguration.hpp"
#include "utils/assertion.hpp"
#include "xml/ConfigParser.hpp"
#include "xml/XMLAttribute.hpp"

namespace precice::action {

ActionConfiguration::ActionConfiguration(
    xml::XMLTag               &parent,
    mesh::PtrMeshConfiguration meshConfig)
    : NAME_DIVIDE_BY_AREA("divide-by-area"),
      NAME_MULTIPLY_BY_AREA("multiply-by-area"),
      NAME_SUMMATION("summation"),
      NAME_PYTHON("python"),
      NAME_RECORDER("recorder"),
      TAG_SOURCE_DATA("source-data"),
      TAG_TARGET_DATA("target-data"),
      TAG_CONVERGENCE_TOLERANCE("convergence-tolerance"),
      TAG_MAX_ITERATIONS("max-iterations"),
      TAG_MODULE_PATH("path"),
      TAG_MODULE_NAME("module"),
      WRITE_MAPPING_POST("write-mapping-post"),
      READ_MAPPING_POST("read-mapping-post"),
      _meshConfig(std::move(meshConfig))
{
  using namespace xml;
  XMLTag tagSourceData(*this, TAG_SOURCE_DATA, XMLTag::OCCUR_ONCE);
  tagSourceData.setDocumentation("Single data to read from. ");
  XMLTag tagMultipleSourceData(*this, TAG_SOURCE_DATA, XMLTag::OCCUR_ONCE_OR_MORE);
  tagMultipleSourceData.setDocumentation("Multiple data to read from.");
  XMLTag tagTargetData(*this, TAG_TARGET_DATA, XMLTag::OCCUR_ONCE);
  tagTargetData.setDocumentation("Data to read from and write to.");

  auto attrName = XMLAttribute<std::string>(ATTR_NAME).setDocumentation("Name of the data.");
  tagSourceData.addAttribute(attrName);
  tagMultipleSourceData.addAttribute(attrName);
  tagTargetData.addAttribute(attrName);

  std::list<XMLTag>  tags;
  XMLTag::Occurrence occ = XMLTag::OCCUR_ARBITRARY;
  {
    XMLTag tag(*this, NAME_MULTIPLY_BY_AREA, occ, TAG);
    tag.setDocumentation("Multiplies data values with mesh area associated to vertex holding the value.");
    tag.addSubtag(tagTargetData);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, NAME_DIVIDE_BY_AREA, occ, TAG);
    tag.setDocumentation("Divides data values by mesh area associated to vertex holding the value.");
    tag.addSubtag(tagTargetData);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, NAME_SUMMATION, occ, TAG);
    tag.setDocumentation("Sums up multiple source data values and writes the result into target data.");
    tag.addSubtag(tagMultipleSourceData);
    tag.addSubtag(tagTargetData);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, NAME_RECORDER, occ, TAG);
    tag.setDocumentation("Records action invocations for testing purposes.");
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, NAME_PYTHON, occ, TAG);
    tag.setDocumentation("Calls Python script to execute action."
                         " See preCICE file \"src/action/PythonAction.py\" for an example.");

    XMLTag tagModulePath(*this, TAG_MODULE_PATH, XMLTag::OCCUR_NOT_OR_ONCE);
    tagModulePath.setDocumentation("Directory path to Python module, i.e. script file."
                                   " If it doesn't occur, the current path is used");
    tagModulePath.addAttribute(makeXMLAttribute(ATTR_NAME, "").setDocumentation("The path to the directory of the module."));
    tag.addSubtag(tagModulePath);

    XMLTag tagModule(*this, TAG_MODULE_NAME, XMLTag::OCCUR_ONCE);
    tagModule.setDocumentation("Name of Python module, i.e. Python script file without file ending. "
                               "The module name has to differ from existing (library) modules, "
                               "otherwise, the existing module will be loaded instead of the user script.");
    tagModule.addAttribute(attrName);
    tag.addSubtag(tagModule);

    XMLTag tagOptionalSourceData(*this, TAG_SOURCE_DATA, XMLTag::OCCUR_NOT_OR_ONCE);
    tagOptionalSourceData.setDocumentation("Source data to be read is handed to the Python module."
                                           " Can be omitted, if only a target data is needed.");
    tagOptionalSourceData.addAttribute(attrName);
    tag.addSubtag(tagOptionalSourceData);

    XMLTag tagOptionalTargetData(*this, TAG_TARGET_DATA, XMLTag::OCCUR_NOT_OR_ONCE);
    tagOptionalTargetData.setDocumentation("Target data to be read and written to is handed to the Python module."
                                           " Can be omitted, if only source data is needed.");
    tagOptionalTargetData.addAttribute(attrName);
    tag.addSubtag(tagOptionalTargetData);

    tags.push_back(tag);
  }

  auto attrTiming = XMLAttribute<std::string>(ATTR_TIMING)
                        .setDocumentation("Determines when (relative to advancing the coupling scheme and the data mappings) the action is executed.")
                        .setOptions({WRITE_MAPPING_POST, READ_MAPPING_POST});

  auto attrMesh = XMLAttribute<std::string>(ATTR_MESH)
                      .setDocumentation("Determines mesh used in action.");
  for (XMLTag &tag : tags) {
    tag.addAttribute(attrTiming);
    tag.addAttribute(attrMesh);
    parent.addSubtag(tag);
  }
}

void ActionConfiguration::xmlTagCallback(
    const xml::ConfigurationContext &context,
    xml::XMLTag                     &callingTag)
{
  PRECICE_TRACE(callingTag.getName());
  if (callingTag.getNamespace() == TAG) {
    _configuredAction        = ConfiguredAction();
    _configuredAction.type   = callingTag.getName();
    _configuredAction.timing = callingTag.getStringAttributeValue(ATTR_TIMING);
    _configuredAction.mesh   = callingTag.getStringAttributeValue(ATTR_MESH);
    // addSubtags ( callingTag, _configured.type );
  } else if (callingTag.getName() == TAG_SOURCE_DATA) {
    _configuredAction.sourceDataVector.push_back(callingTag.getStringAttributeValue(ATTR_NAME));
  } else if (callingTag.getName() == TAG_TARGET_DATA) {
    _configuredAction.targetData = callingTag.getStringAttributeValue(ATTR_NAME);
  } else if (callingTag.getName() == TAG_CONVERGENCE_TOLERANCE) {
    _configuredAction.convergenceTolerance =
        callingTag.getDoubleAttributeValue(ATTR_VALUE);
  } else if (callingTag.getName() == TAG_MAX_ITERATIONS) {
    _configuredAction.maxIterations = callingTag.getIntAttributeValue(ATTR_VALUE);
  } else if (callingTag.getName() == TAG_MODULE_PATH) {
    _configuredAction.path = callingTag.getStringAttributeValue(ATTR_NAME);
  } else if (callingTag.getName() == TAG_MODULE_NAME) {
    _configuredAction.module = callingTag.getStringAttributeValue(ATTR_NAME);
  }
}

void ActionConfiguration::xmlEndTagCallback(
    const xml::ConfigurationContext &context,
    xml::XMLTag                     &callingTag)
{
  if (callingTag.getNamespace() == TAG) {
    createAction();
  }
}

int ActionConfiguration::getUsedMeshID() const
{
  PRECICE_CHECK(_meshConfig->hasMeshName(_configuredAction.mesh), "No mesh name \"{}\" found. Please check that the correct mesh name is used.", _configuredAction.mesh);
  return _meshConfig->getMesh(_configuredAction.mesh)->getID();
}

void ActionConfiguration::createAction()
{
  PRECICE_TRACE();

  PRECICE_ASSERT(_configuredAction.type != std::string(""));
  action::Action::Timing timing = getTiming();

  // Determine data and mesh
  std::vector<int> sourceDataIDs;
  int              targetDataID = -1;
  PRECICE_CHECK(_meshConfig->hasMeshName(_configuredAction.mesh),
                "Data action uses mesh \"{}\" which is not configured. Please ensure that the correct mesh name is given in <action:python mesh=\"...\">", _configuredAction.mesh);
  mesh::PtrMesh mesh = _meshConfig->getMesh(_configuredAction.mesh);

  if (!_configuredAction.targetData.empty()) {
    PRECICE_CHECK(mesh->hasDataName(_configuredAction.targetData),
                  "Data action uses target data \"{}\" which is not configured. Please ensure that the target data name is used by the mesh with name \"{}\".", _configuredAction.targetData, _configuredAction.mesh);
    targetDataID = mesh->data(_configuredAction.targetData)->getID();
    PRECICE_ASSERT(targetDataID != -1);
  }

  for (const std::string &dataName : _configuredAction.sourceDataVector) {
    PRECICE_CHECK(mesh->hasDataName(dataName), "Data action uses source data \"{}\" which is not configured. Please ensure that the target data name is used by the mesh with name \"{}\".", dataName, _configuredAction.mesh);
    sourceDataIDs.push_back(mesh->data(dataName)->getID());
  }

  PRECICE_CHECK((_configuredAction.sourceDataVector.empty() || not sourceDataIDs.empty()),
                "Data action uses source data \"{}\" which is not configured. Please ensure that the source data name is used by the mesh with name \"{}\".", _configuredAction.sourceDataVector.back(), _configuredAction.mesh);

  action::PtrAction action;
  if (_configuredAction.type == NAME_MULTIPLY_BY_AREA) {
    action = action::PtrAction(
        new action::ScaleByAreaAction(timing, targetDataID,
                                      mesh, action::ScaleByAreaAction::SCALING_MULTIPLY_BY_AREA));
  } else if (_configuredAction.type == NAME_DIVIDE_BY_AREA) {
    action = action::PtrAction(
        new action::ScaleByAreaAction(timing, targetDataID,
                                      mesh, action::ScaleByAreaAction::SCALING_DIVIDE_BY_AREA));
  } else if (_configuredAction.type == NAME_SUMMATION) {
    action = action::PtrAction(
        new action::SummationAction(timing, sourceDataIDs, targetDataID, mesh));
  } else if (_configuredAction.type == NAME_RECORDER) {
    action = action::PtrAction(
        new action::RecorderAction(timing, mesh));
  }
#ifndef PRECICE_NO_PYTHON
  else if (_configuredAction.type == NAME_PYTHON) {
    action = action::PtrAction(
        new action::PythonAction(timing, _configuredAction.path, _configuredAction.module,
                                 mesh, targetDataID, sourceDataIDs.back()));
  }
#endif
  PRECICE_ASSERT(action.get() != nullptr);
  _actions.push_back(std::move(action));
}

action::Action::Timing ActionConfiguration::getTiming() const
{
  PRECICE_TRACE(_configuredAction.timing);
  action::Action::Timing timing;
  if (_configuredAction.timing == WRITE_MAPPING_POST) {
    timing = action::Action::WRITE_MAPPING_POST;
  } else if (_configuredAction.timing == READ_MAPPING_POST) {
    timing = action::Action::READ_MAPPING_POST;
  } else {
    PRECICE_ERROR("Unknown action timing \"{}\". "
                  "Valid action timings are read-mapping-post and write-mapping-post.",
                  _configuredAction.timing);
  }
  return timing;
}

} // namespace precice::action
