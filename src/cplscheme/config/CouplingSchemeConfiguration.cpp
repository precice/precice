// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "CouplingSchemeConfiguration.hpp"
#include "cplscheme/config/PostProcessingConfiguration.hpp"
#include "cplscheme/SerialCouplingScheme.hpp"
#include "cplscheme/ParallelCouplingScheme.hpp"
#include "cplscheme/MultiCouplingScheme.hpp"
#include "cplscheme/CompositionalCouplingScheme.hpp"
#include "cplscheme/impl/ConvergenceMeasure.hpp"
#include "cplscheme/impl/AbsoluteConvergenceMeasure.hpp"
#include "cplscheme/impl/RelativeConvergenceMeasure.hpp"
#include "cplscheme/impl/ResidualRelativeConvergenceMeasure.hpp"
#include "cplscheme/impl/MinIterationConvergenceMeasure.hpp"
#include "cplscheme/impl/PostProcessing.hpp"
#include "cplscheme/SharedPointer.hpp"
#include "mesh/config/MeshConfiguration.hpp"
#include "mesh/config/DataConfiguration.hpp"
#include "m2n/M2N.hpp"
#include "m2n/config/M2NConfiguration.hpp"
#include "utils/Globals.hpp"
#include "utils/Dimensions.hpp"
#include "utils/xml/XMLAttribute.hpp"
#include "utils/xml/ValidatorEquals.hpp"
#include "utils/xml/ValidatorOr.hpp"
#include "utils/xml/XMLTag.hpp"
#include "precice/impl/SharedPointer.hpp"
#include "precice/impl/Participant.hpp"

namespace precice {
namespace cplscheme {

using precice::impl::PtrParticipant;

tarch::logging::Log CouplingSchemeConfiguration::
   _log("precice::cplscheme::CouplingSchemeConfiguration");

CouplingSchemeConfiguration:: CouplingSchemeConfiguration
(
  utils::XMLTag&                            parent,
  const mesh::PtrMeshConfiguration&         meshConfig,
  const m2n::M2NConfiguration::SharedPointer& m2nConfig)
:
  TAG("coupling-scheme"),
  TAG_PARTICIPANTS("participants"),
  TAG_PARTICIPANT("participant"),
  TAG_EXCHANGE("exchange"),
  TAG_MAX_TIME("max-time"),
  TAG_MAX_TIMESTEPS("max-timesteps"),
  TAG_TIMESTEP_LENGTH("timestep-length"),
  TAG_ABS_CONV_MEASURE("absolute-convergence-measure"),
  TAG_REL_CONV_MEASURE("relative-convergence-measure"),
  TAG_RES_REL_CONV_MEASURE("residual-relative-convergence-measure"),
  TAG_MIN_ITER_CONV_MEASURE("min-iteration-convergence-measure"),
  TAG_MAX_ITERATIONS("max-iterations"),
  TAG_CHECKPOINT("checkpoint"),
  TAG_EXTRAPOLATION("extrapolation-order"),
  ATTR_DATA("data"),
  ATTR_MESH("mesh"),
  ATTR_PARTICIPANT("participant"),
  ATTR_INITIALIZE("initialize"),
  ATTR_TYPE("type"),
  ATTR_FIRST("first"),
  ATTR_SECOND("second"),
  ATTR_VALUE("value"),
  ATTR_VALID_DIGITS("valid-digits"),
  ATTR_METHOD("method"),
  ATTR_LIMIT("limit"),
  ATTR_MIN_ITERATIONS("min-iterations"),
  ATTR_NAME("name"),
  ATTR_TIMESTEP_INTERVAL("timestep-interval"),
  ATTR_FROM("from"),
  ATTR_TO("to"),
  ATTR_SUFFICES("suffices"),
  ATTR_CONTROL("control"),
  VALUE_SERIAL_EXPLICIT("serial-explicit"),
  VALUE_PARALLEL_EXPLICIT("parallel-explicit"),
  VALUE_SERIAL_IMPLICIT("serial-implicit"),
  VALUE_PARALLEL_IMPLICIT("parallel-implicit"),
  VALUE_MULTI("multi"),
  VALUE_UNCOUPLED("uncoupled"),
  VALUE_FIXED("fixed"),
  VALUE_FIRST_PARTICIPANT("first-participant"),
  _config(),
  _meshConfig(meshConfig),
  _m2nConfig(m2nConfig),
  _isValid(false),
  _couplingSchemes(),
  _couplingSchemeCompositions()
{
  using namespace utils;
  //preciceCheck ( _couplingSchemes.size() == 0, "parseSubtag()",
  //               "Only one tag <coupling-scheme> can be defined!" );

  XMLTag::Occurrence occ = XMLTag::OCCUR_ARBITRARY;
  std::list<XMLTag> tags;
  std::string doc;

  {
    XMLTag tag(*this, VALUE_SERIAL_EXPLICIT, occ, TAG);
    doc = "Explicit coupling scheme according to conventional serial";
    doc += " staggered procedure (CSS).";
    tag.setDocumentation(doc);
    addTypespecifcSubtags(VALUE_SERIAL_EXPLICIT, tag);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_PARALLEL_EXPLICIT, occ, TAG);
    doc = "Explicit coupling scheme according to conventional parallel";
    doc += " staggered procedure (CPS).";
    tag.setDocumentation(doc);
    addTypespecifcSubtags(VALUE_PARALLEL_EXPLICIT, tag);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_SERIAL_IMPLICIT, occ, TAG);
    doc = "Implicit coupling scheme according to block Gauss-Seidel iterations (S-System).";
    doc += " Improved implicit iterations are achieved by using a post-processing.";
    tag.setDocumentation(doc);
    addTypespecifcSubtags(VALUE_SERIAL_IMPLICIT, tag);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_PARALLEL_IMPLICIT, occ, TAG);
    doc = "Parallel Implicit coupling scheme according to block Jacobi iterations (V-System).";
    doc += " Improved implicit iterations are achieved by using a post-processing.";
    tag.setDocumentation(doc);
    addTypespecifcSubtags(VALUE_PARALLEL_IMPLICIT, tag);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_MULTI, occ, TAG);
    doc = "Multi coupling scheme according to block Jacobi iterations.";
    doc += " Improved implicit iterations are achieved by using a post-processing (recommended!).";
    tag.setDocumentation(doc);
    addTypespecifcSubtags(VALUE_MULTI, tag);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_UNCOUPLED, occ, TAG);
    doc = "Coupling scheme for using geometry mode, i.e., for a solver that uses";
    doc += " preCICE without coupling to another solver.";
    tag.setDocumentation(doc);
    addTypespecifcSubtags(VALUE_UNCOUPLED, tag);
    tags.push_back(tag);
  }

  foreach (XMLTag& tag, tags){
    parent.addSubtag(tag);
  }
}

bool CouplingSchemeConfiguration:: hasCouplingScheme
(
  const std::string& participantName ) const
{
  return utils::contained(participantName, _couplingSchemes);
}

const PtrCouplingScheme& CouplingSchemeConfiguration:: getCouplingScheme
(
  const std::string& participantName ) const
{
  preciceCheck(utils::contained(participantName, _couplingSchemes),
               "getCouplingScheme()", "No coupling scheme defined for "
               << "participant \"" << participantName << "\"!");
  return _couplingSchemes.find(participantName)->second;
}

void CouplingSchemeConfiguration:: xmlTagCallback
(
  utils::XMLTag& tag )
{
  preciceTrace1("xmlTagCallback()", tag.getFullName());
  if (tag.getNamespace() == TAG){
    _config.type = tag.getName();
    _postProcConfig->clear();
  }
  else if (tag.getName() == TAG_CHECKPOINT){
    _config.checkpointTimestepInterval =
        tag.getIntAttributeValue(ATTR_TIMESTEP_INTERVAL);
  }
  else if (tag.getName() == TAG_PARTICIPANTS){
    assertion(_config.type != VALUE_UNCOUPLED);
    _config.participants.push_back(tag.getStringAttributeValue(ATTR_FIRST));
    _config.participants.push_back(tag.getStringAttributeValue(ATTR_SECOND));
  }
  else if (tag.getName() == TAG_PARTICIPANT){
    assertion(_config.type == VALUE_MULTI);
    bool control = tag.getBooleanAttributeValue(ATTR_CONTROL);
    if(control){
      preciceCheck ( not _config.setController, "xmlTagCallback()",
                       "Only one controller per MultiCoupling can be defined" );
      _config.controller = tag.getStringAttributeValue(ATTR_NAME);
      _config.setController = true;
    }
    else{
      _config.participants.push_back(tag.getStringAttributeValue(ATTR_NAME));
    }

  }
  else if (tag.getName() == TAG_MAX_TIME){
    _config.maxTime = tag.getDoubleAttributeValue(ATTR_VALUE);
  }
  else if (tag.getName() == TAG_MAX_TIMESTEPS){
    _config.maxTimesteps =
        tag.getIntAttributeValue(ATTR_VALUE);
  }
  else if ( tag.getName() == TAG_TIMESTEP_LENGTH ) {
    _config.timestepLength =
        tag.getDoubleAttributeValue(ATTR_VALUE);
    _config.validDigits =
        tag.getIntAttributeValue(ATTR_VALID_DIGITS);
    _config.dtMethod = getTimesteppingMethod(
        tag.getStringAttributeValue(ATTR_METHOD));
  }
  else if ( tag.getName() == TAG_ABS_CONV_MEASURE ) {
    std::string dataName = tag.getStringAttributeValue(ATTR_DATA);
    std::string meshName = tag.getStringAttributeValue(ATTR_MESH);
    double limit = tag.getDoubleAttributeValue(ATTR_LIMIT);
    bool suffices = tag.getBooleanAttributeValue(ATTR_SUFFICES);
    assertion(_config.type == VALUE_SERIAL_IMPLICIT || _config.type == VALUE_PARALLEL_IMPLICIT
        || _config.type == VALUE_MULTI);
    addAbsoluteConvergenceMeasure(dataName, meshName, limit, suffices);
  }
  else if ( tag.getName() == TAG_REL_CONV_MEASURE ) {
    std::string dataName = tag.getStringAttributeValue (ATTR_DATA);
    std::string meshName = tag.getStringAttributeValue(ATTR_MESH);
    double limit = tag.getDoubleAttributeValue(ATTR_LIMIT);
    bool suffices = tag.getBooleanAttributeValue(ATTR_SUFFICES);
    assertion(_config.type == VALUE_SERIAL_IMPLICIT || _config.type == VALUE_PARALLEL_IMPLICIT
        || _config.type == VALUE_MULTI);
    addRelativeConvergenceMeasure(dataName, meshName, limit, suffices);
  }
  else if ( tag.getName() == TAG_RES_REL_CONV_MEASURE ) {
    std::string dataName = tag.getStringAttributeValue (ATTR_DATA);
    std::string meshName = tag.getStringAttributeValue(ATTR_MESH);
    double limit = tag.getDoubleAttributeValue(ATTR_LIMIT);
    bool suffices = tag.getBooleanAttributeValue(ATTR_SUFFICES);
    assertion(_config.type == VALUE_SERIAL_IMPLICIT || _config.type == VALUE_PARALLEL_IMPLICIT
        || _config.type == VALUE_MULTI);
    addResidualRelativeConvergenceMeasure(dataName, meshName, limit, suffices);
  }
  else if ( tag.getName() == TAG_MIN_ITER_CONV_MEASURE ) {
    std::string dataName = tag.getStringAttributeValue(ATTR_DATA);
    std::string meshName = tag.getStringAttributeValue(ATTR_MESH);
    int minIterations = tag.getIntAttributeValue(ATTR_MIN_ITERATIONS);
    bool suffices = tag.getBooleanAttributeValue(ATTR_SUFFICES);
    assertion(_config.type == VALUE_SERIAL_IMPLICIT || _config.type == VALUE_PARALLEL_IMPLICIT
        || _config.type == VALUE_MULTI);
    addMinIterationConvergenceMeasure(dataName, meshName, minIterations, suffices);
  }
  else if (tag.getName() == TAG_EXCHANGE){
    assertion(_config.type != VALUE_UNCOUPLED);
    std::string nameData = tag.getStringAttributeValue(ATTR_DATA);
    std::string nameMesh = tag.getStringAttributeValue(ATTR_MESH);
    std::string nameParticipantFrom = tag.getStringAttributeValue(ATTR_FROM);
    std::string nameParticipantTo = tag.getStringAttributeValue(ATTR_TO);
    bool initialize = tag.getBooleanAttributeValue(ATTR_INITIALIZE);
    mesh::PtrData exchangeData;
    mesh::PtrMesh exchangeMesh;
    foreach (mesh::PtrMesh mesh, _meshConfig->meshes()){
      if ( mesh->getName() == nameMesh ) {
        for ( mesh::PtrData data : mesh->data() ) {
          if ( data->getName() == nameData ) {
            exchangeData = data;
            exchangeMesh = mesh;
            break;
          }
        }
      }
    }
    if (exchangeData.get() == NULL){
      std::ostringstream stream;
      stream << "Mesh \"" << nameMesh << "\" with data \"" << nameData
             << "\" not defined at definition of coupling scheme";
      throw stream.str();
    }
    _meshConfig->addNeededMesh(nameParticipantFrom, nameMesh);
    _meshConfig->addNeededMesh(nameParticipantTo, nameMesh);
    _config.exchanges.push_back(boost::make_tuple(exchangeData, exchangeMesh,
                  nameParticipantFrom,nameParticipantTo, initialize));
  }
  else if (tag.getName() == TAG_MAX_ITERATIONS){
    assertion(_config.type == VALUE_SERIAL_IMPLICIT || _config.type == VALUE_PARALLEL_IMPLICIT
        || _config.type == VALUE_MULTI);
    _config.maxIterations = tag.getIntAttributeValue(ATTR_VALUE);
  }
  else if (tag.getName() == TAG_EXTRAPOLATION){
    assertion(_config.type == VALUE_SERIAL_IMPLICIT || _config.type == VALUE_PARALLEL_IMPLICIT
        || _config.type == VALUE_MULTI);
    _config.extrapolationOrder = tag.getIntAttributeValue(ATTR_VALUE);
  }
}

void CouplingSchemeConfiguration:: xmlEndTagCallback
(
  utils::XMLTag& tag )
{
  preciceTrace1("xmlEndTagCallback()", tag.getFullName());
  if (tag.getNamespace() == TAG){
    if (_config.type == VALUE_SERIAL_EXPLICIT){
      std::string accessor(_config.participants[0]);
      PtrCouplingScheme scheme = createSerialExplicitCouplingScheme(accessor);
      addCouplingScheme(scheme, accessor);
      //_couplingSchemes[accessor] = scheme;
      accessor = _config.participants[1];
      scheme = createSerialExplicitCouplingScheme(accessor);
      addCouplingScheme(scheme, accessor);
      //_couplingSchemes[accessor] = scheme;
      _config = Config();
    }
    else if (_config.type == VALUE_PARALLEL_EXPLICIT){
      std::string accessor(_config.participants[0]);
      PtrCouplingScheme scheme = createParallelExplicitCouplingScheme(accessor);
      addCouplingScheme(scheme, accessor);
      //_couplingSchemes[accessor] = scheme;
      accessor = _config.participants[1];
      scheme = createParallelExplicitCouplingScheme(accessor);
      addCouplingScheme(scheme, accessor);
      //_couplingSchemes[accessor] = scheme;
      _config = Config();
    }
    else if (_config.type == VALUE_SERIAL_IMPLICIT){
      std::string accessor(_config.participants[0]);
      PtrCouplingScheme scheme = createSerialImplicitCouplingScheme(accessor);
      addCouplingScheme(scheme, accessor);
      //_couplingSchemes[accessor] = scheme;
      accessor = _config.participants[1];
      scheme = createSerialImplicitCouplingScheme(accessor);
      addCouplingScheme(scheme, accessor);
      //_couplingSchemes[accessor] = scheme;
      _config = Config();
    }
    else if (_config.type == VALUE_PARALLEL_IMPLICIT){
      std::string accessor(_config.participants[0]);
      PtrCouplingScheme scheme = createParallelImplicitCouplingScheme(accessor);
      addCouplingScheme(scheme, accessor);
      accessor = _config.participants[1];
      scheme = createParallelImplicitCouplingScheme(accessor);
      addCouplingScheme(scheme, accessor);
      _config = Config();
    }
    else if (_config.type == VALUE_MULTI){
      preciceCheck ( _config.setController, "xmlTagCallback()",
              "One controller per MultiCoupling need to be defined" );
      foreach(std::string& accessor, _config.participants){
        PtrCouplingScheme scheme = createMultiCouplingScheme(accessor);
        addCouplingScheme(scheme, accessor);
      }
      PtrCouplingScheme scheme = createMultiCouplingScheme(_config.controller);
      addCouplingScheme(scheme, _config.controller);
      _config = Config();
    }
    else if (_config.type == VALUE_UNCOUPLED){
      assertion(false);
    }
    else {
      assertion1(false,_config.type);
    }
  }
}

void CouplingSchemeConfiguration:: addCouplingScheme
(
  PtrCouplingScheme  cplScheme,
  const std::string& participantName )
{
  preciceTrace1 ( "addCouplingScheme()", participantName );
  if (utils::contained(participantName, _couplingSchemes)) {
    preciceDebug("Coupling scheme exists already for participant");
    if (utils::contained(participantName, _couplingSchemeCompositions)) {
      preciceDebug("Coupling scheme composition exists already for participant");
      // Fetch the composition and add the new scheme.
      assertion(_couplingSchemeCompositions[participantName] != NULL);
      _couplingSchemeCompositions[participantName]->addCouplingScheme(cplScheme);
    }
    else {
      preciceDebug("No composition exists for the participant");
      // No composition exists, thus, the existing scheme is no composition.
      // Create a new composition, add the already existing and new scheme, and
      // overwrite the existing scheme with the composition.
      CompositionalCouplingScheme* composition = new CompositionalCouplingScheme();
      composition->addCouplingScheme(_couplingSchemes[participantName]);
      composition->addCouplingScheme(cplScheme);
      _couplingSchemes[participantName] = PtrCouplingScheme(composition);
    }
  }
  else {
    preciceDebug("No coupling scheme exists for the participant");
    // Store the new coupling scheme.
    _couplingSchemes[participantName] = cplScheme;
  }
}

void CouplingSchemeConfiguration:: addTypespecifcSubtags
(
  const std::string& type,
  //const std::string& name,
  utils::XMLTag&     tag  )
{
  preciceTrace1( "addTypespecifcSubtags()", type );
  addTransientLimitTags(tag);
  _config.type = type;
  //_config.name = name;

  addTagCheckpoint(tag);

  if (type == VALUE_SERIAL_EXPLICIT){
    addTagParticipants(tag);
    addTagExchange(tag);
  }
  else if (type == VALUE_PARALLEL_EXPLICIT){
    addTagParticipants(tag);
    addTagExchange(tag);
  }
  else if ( type == VALUE_PARALLEL_IMPLICIT ) {
    addTagParticipants(tag);
    addTagExchange(tag);
    addTagAbsoluteConvergenceMeasure(tag);
    addTagRelativeConvergenceMeasure(tag);
    addTagResidualRelativeConvergenceMeasure(tag);
    addTagMinIterationConvergenceMeasure(tag);
    addTagMaxIterations(tag);
    addTagExtrapolation(tag);
    addTagPostProcessing(tag);
  }
  else if ( type == VALUE_MULTI ) {
    addTagParticipant(tag);
    addTagExchange(tag);
    addTagAbsoluteConvergenceMeasure(tag);
    addTagRelativeConvergenceMeasure(tag);
    addTagResidualRelativeConvergenceMeasure(tag);
    addTagMinIterationConvergenceMeasure(tag);
    addTagMaxIterations(tag);
    addTagExtrapolation(tag);
    addTagPostProcessing(tag);
  }
  else if ( type == VALUE_SERIAL_IMPLICIT ) {
    addTagParticipants(tag);
    addTagExchange(tag);
    addTagAbsoluteConvergenceMeasure(tag);
    addTagRelativeConvergenceMeasure(tag);
    addTagResidualRelativeConvergenceMeasure(tag);
    addTagMinIterationConvergenceMeasure(tag);
    addTagMaxIterations(tag);
    addTagExtrapolation(tag);
    addTagPostProcessing(tag);
  }
  else if (type == VALUE_UNCOUPLED){
  }
  else {
    preciceError("addTypespecificSubtags()", "Unknown coupling scheme type!");
  }
}

void CouplingSchemeConfiguration:: addTagCheckpoint
(
   utils::XMLTag& tag )
{
  using namespace utils;
  XMLTag tagCheckpoint(*this, TAG_CHECKPOINT, XMLTag::OCCUR_NOT_OR_ONCE);
   utils::XMLAttribute<int> attrTimestepInterval(ATTR_TIMESTEP_INTERVAL);
//   utils::ValidatorGreaterThan<int> validTimestepInterval ( 0 );
//   attrTimestepInterval.setValidator ( validTimestepInterval );
   tagCheckpoint.addAttribute(attrTimestepInterval);
   tag.addSubtag(tagCheckpoint);
}

void CouplingSchemeConfiguration:: addTransientLimitTags
(
   utils::XMLTag & tag )
{
  using namespace utils;
  XMLTag tagMaxTime(*this, TAG_MAX_TIME, XMLTag::OCCUR_NOT_OR_ONCE);
  XMLAttribute<double> attrValueMaxTime(ATTR_VALUE);
  tagMaxTime.addAttribute(attrValueMaxTime);
  tag.addSubtag(tagMaxTime);

  XMLTag tagMaxTimesteps(*this, TAG_MAX_TIMESTEPS, XMLTag::OCCUR_NOT_OR_ONCE);
  XMLAttribute<int> attrValueMaxTimesteps(ATTR_VALUE);
  tagMaxTimesteps.addAttribute(attrValueMaxTimesteps);
  tag.addSubtag(tagMaxTimesteps);

  XMLTag tagTimestepLength(*this, TAG_TIMESTEP_LENGTH, XMLTag::OCCUR_ONCE );
  XMLAttribute<double> attrValueTimestepLength(ATTR_VALUE);
  attrValueTimestepLength.setDefaultValue(CouplingScheme::UNDEFINED_TIMESTEP_LENGTH);
  tagTimestepLength.addAttribute(attrValueTimestepLength);
  XMLAttribute<int> attrValidDigits(ATTR_VALID_DIGITS);
  attrValidDigits.setDefaultValue(10);
  tagTimestepLength.addAttribute(attrValidDigits);
  XMLAttribute<std::string> attrMethod(ATTR_METHOD);
  attrMethod.setDefaultValue(VALUE_FIXED);
  ValidatorEquals<std::string> validFixed(VALUE_FIXED);
  ValidatorEquals<std::string> validFirst(VALUE_FIRST_PARTICIPANT);
//  ValidatorEquals<std::string> validSec ( TagTimestepLength::VALUE_SECOND_PARTICIPANT );
  attrMethod.setValidator(validFixed || validFirst);
  tagTimestepLength.addAttribute(attrMethod);
  tag.addSubtag(tagTimestepLength);
}

void CouplingSchemeConfiguration:: addTagParticipants
(
  utils::XMLTag& tag )
{
  using namespace utils;
  XMLTag tagParticipants(*this, TAG_PARTICIPANTS, XMLTag::OCCUR_ONCE);
  XMLAttribute<std::string> attrFirst(ATTR_FIRST);
  tagParticipants.addAttribute(attrFirst);
  XMLAttribute<std::string> attrSecond(ATTR_SECOND);
  tagParticipants.addAttribute(attrSecond);
  tag.addSubtag(tagParticipants);
}

void CouplingSchemeConfiguration:: addTagParticipant
(
  utils::XMLTag& tag )
{
  using namespace utils;
  XMLTag tagParticipant(*this, TAG_PARTICIPANT, XMLTag::OCCUR_ONCE_OR_MORE);
  XMLAttribute<std::string> attrName(ATTR_NAME);
  tagParticipant.addAttribute(attrName);
  XMLAttribute<bool> attrControl(ATTR_CONTROL);
  attrControl.setDefaultValue(false);
  tagParticipant.addAttribute(attrControl);
  tag.addSubtag(tagParticipant);
}

void CouplingSchemeConfiguration:: addTagExchange
(
  utils::XMLTag& tag )
{
  using namespace utils;
  XMLTag tagExchange(*this, TAG_EXCHANGE, XMLTag::OCCUR_ONCE_OR_MORE);
  XMLAttribute<std::string> attrData(ATTR_DATA);
  tagExchange.addAttribute(attrData);
  XMLAttribute<std::string> attrMesh(ATTR_MESH);
  tagExchange.addAttribute(attrMesh);
  XMLAttribute<std::string> participantFrom(ATTR_FROM);
  tagExchange.addAttribute(participantFrom);
  XMLAttribute<std::string> participantTo(ATTR_TO);
  tagExchange.addAttribute(participantTo);
  XMLAttribute<bool> attrInitialize(ATTR_INITIALIZE);
  attrInitialize.setDefaultValue(false);
  tagExchange.addAttribute(attrInitialize);
  tag.addSubtag(tagExchange);
}

void CouplingSchemeConfiguration:: addTagAbsoluteConvergenceMeasure
(
  utils::XMLTag& tag )
{
  using namespace utils;
  XMLTag tagConvergenceMeasure(*this, TAG_ABS_CONV_MEASURE, XMLTag::OCCUR_ARBITRARY);
  addBaseAttributesTagConvergenceMeasure(tagConvergenceMeasure);
  XMLAttribute<double> attrLimit(ATTR_LIMIT);
  tagConvergenceMeasure.addAttribute(attrLimit);
  tag.addSubtag(tagConvergenceMeasure);
}

void CouplingSchemeConfiguration:: addTagResidualRelativeConvergenceMeasure
(
  utils::XMLTag& tag )
{
  using namespace utils;
  XMLTag tagConvergenceMeasure(*this, TAG_RES_REL_CONV_MEASURE,
                               XMLTag::OCCUR_ARBITRARY );
  addBaseAttributesTagConvergenceMeasure(tagConvergenceMeasure);
  XMLAttribute<double> attrLimit(ATTR_LIMIT);
  tagConvergenceMeasure.addAttribute(attrLimit);
  tag.addSubtag(tagConvergenceMeasure);
}

void CouplingSchemeConfiguration:: addTagRelativeConvergenceMeasure
(
  utils::XMLTag& tag )
{
  using namespace utils;
  XMLTag tagConvergenceMeasure(*this, TAG_REL_CONV_MEASURE, XMLTag::OCCUR_ARBITRARY);
  addBaseAttributesTagConvergenceMeasure(tagConvergenceMeasure);
  XMLAttribute<double> attrLimit(ATTR_LIMIT);
  tagConvergenceMeasure.addAttribute(attrLimit);
  tag.addSubtag(tagConvergenceMeasure);
}

void CouplingSchemeConfiguration:: addTagMinIterationConvergenceMeasure
(
  utils::XMLTag& tag )
{
  utils::XMLTag tagMinIterationConvMeasure (*this,
    TAG_MIN_ITER_CONV_MEASURE, utils::XMLTag::OCCUR_ARBITRARY );
  addBaseAttributesTagConvergenceMeasure(tagMinIterationConvMeasure);
  utils::XMLAttribute<int> attrMinIterations(ATTR_MIN_ITERATIONS);
  tagMinIterationConvMeasure.addAttribute(attrMinIterations);
  tag.addSubtag(tagMinIterationConvMeasure);
}

void CouplingSchemeConfiguration:: addBaseAttributesTagConvergenceMeasure
(
  utils::XMLTag& tag )
{
  utils::XMLAttribute<std::string> attrData(ATTR_DATA);
  attrData.setDocumentation("Data to be measured.");
  tag.addAttribute(attrData);
  utils::XMLAttribute<std::string> attrMesh(ATTR_MESH);
  attrMesh.setDocumentation("Mesh holding the data.");
  tag.addAttribute(attrMesh);
  utils::XMLAttribute<bool> attrSuffices(ATTR_SUFFICES);
  attrSuffices.setDocumentation("If true, suffices to lead to convergence.");
  attrSuffices.setDefaultValue(false);
  tag.addAttribute(attrSuffices);
}

void CouplingSchemeConfiguration:: addTagMaxIterations
(
  utils::XMLTag& tag )
{
  using namespace utils;
  XMLTag tagMaxIterations(*this, TAG_MAX_ITERATIONS, XMLTag::OCCUR_NOT_OR_ONCE);
  XMLAttribute<int> attrValue(ATTR_VALUE);
  tagMaxIterations.addAttribute(attrValue);
  tag.addSubtag(tagMaxIterations);
}

void CouplingSchemeConfiguration:: addTagExtrapolation
(
  utils::XMLTag& tag )
{
  using namespace utils;
  XMLTag tagExtrapolation(*this, TAG_EXTRAPOLATION, XMLTag::OCCUR_NOT_OR_ONCE);
  XMLAttribute<int> attrValue(ATTR_VALUE);
  tagExtrapolation.addAttribute(attrValue);
  tag.addSubtag(tagExtrapolation);
}

void CouplingSchemeConfiguration:: addTagPostProcessing
(
  utils::XMLTag& tag )
{
  preciceTrace1( "addTagPostProcessing()",tag.getFullName());
  if(_postProcConfig.get()==NULL){
    _postProcConfig = PtrPostProcessingConfiguration(
                          new PostProcessingConfiguration(_meshConfig));
  }
  _postProcConfig->connectTags(tag);
}


void CouplingSchemeConfiguration:: addAbsoluteConvergenceMeasure
(
  const std::string& dataName,
  const std::string& meshName,
  double             limit,
  bool               suffices )
{
  preciceTrace( "addAbsoluteConvergenceMeasure()");
  impl::PtrConvergenceMeasure measure(new impl::AbsoluteConvergenceMeasure(limit));
  int dataID = getData(dataName, meshName)->getID();
  _config.convMeasures.push_back(boost::make_tuple(dataID, suffices, meshName, measure));
}

void CouplingSchemeConfiguration:: addRelativeConvergenceMeasure
(
  const std::string& dataName,
  const std::string& meshName,
  double             limit,
  bool               suffices )
{
  preciceTrace( "addRelativeConvergenceMeasure()");
  impl::PtrConvergenceMeasure measure (
      new impl::RelativeConvergenceMeasure(limit) );
  int dataID = getData(dataName, meshName)->getID();
  _config.convMeasures.push_back(boost::make_tuple(dataID, suffices, meshName, measure));
}

void CouplingSchemeConfiguration:: addResidualRelativeConvergenceMeasure
(
  const std::string& dataName,
  const std::string& meshName,
  double             limit,
  bool               suffices )
{
  preciceTrace( "addResidualRelativeConvergenceMeasure()");
  impl::PtrConvergenceMeasure measure (
      new impl::ResidualRelativeConvergenceMeasure(limit) );
  int dataID = getData(dataName, meshName)->getID();
  _config.convMeasures.push_back(boost::make_tuple(dataID, suffices, meshName, measure));
}

void CouplingSchemeConfiguration:: addMinIterationConvergenceMeasure
(
  const std::string& dataName,
  const std::string& meshName,
  int                minIterations,
  bool               suffices )
{
  preciceTrace( "addMinIterationConvergenceMeasure()");
  impl::PtrConvergenceMeasure measure (
      new impl::MinIterationConvergenceMeasure(minIterations) );
  int dataID = getData(dataName, meshName)->getID();
  _config.convMeasures.push_back(boost::make_tuple(dataID, suffices, meshName, measure));
}

mesh::PtrData CouplingSchemeConfiguration:: getData
(
  const std::string& dataName,
  const std::string& meshName ) const
{
  foreach ( mesh::PtrMesh mesh, _meshConfig->meshes() ){
    if ( meshName == mesh->getName() ){
      for ( mesh::PtrData data : mesh->data() ){
        if ( dataName == data->getName() ){
          return data;
        }
      }
    }
  }
  preciceError ( "getData()", "Data \"" << dataName << "\" used by mesh \""
                 << meshName << "\" is not configured!" );
}


PtrCouplingScheme CouplingSchemeConfiguration:: createSerialExplicitCouplingScheme
(
  const std::string& accessor ) const
{
  preciceTrace1("createSerialExplicitCouplingScheme()", accessor);
  //assertion ( not utils::contained(accessor, _couplingSchemes) );
  m2n::M2N::SharedPointer m2n = _m2nConfig->getM2N (
      _config.participants[0], _config.participants[1] );
  SerialCouplingScheme* scheme = new SerialCouplingScheme (
      _config.maxTime, _config.maxTimesteps, _config.timestepLength,
      _config.validDigits, _config.participants[0], _config.participants[1],
      accessor, m2n, _config.dtMethod, BaseCouplingScheme::Explicit );
  scheme->setCheckPointTimestepInterval ( _config.checkpointTimestepInterval );

  addDataToBeExchanged(*scheme, accessor);

  return PtrCouplingScheme(scheme);
}

PtrCouplingScheme CouplingSchemeConfiguration:: createParallelExplicitCouplingScheme
(
  const std::string& accessor ) const
{
  preciceTrace1("createParallelExplicitCouplingScheme()", accessor);
  //assertion ( not utils::contained(accessor, _couplingSchemes) );
  m2n::M2N::SharedPointer m2n = _m2nConfig->getM2N (
      _config.participants[0], _config.participants[1] );
  ParallelCouplingScheme* scheme = new ParallelCouplingScheme (
      _config.maxTime, _config.maxTimesteps, _config.timestepLength,
      _config.validDigits, _config.participants[0], _config.participants[1],
      accessor, m2n, _config.dtMethod, BaseCouplingScheme::Explicit );
  scheme->setCheckPointTimestepInterval ( _config.checkpointTimestepInterval );

  addDataToBeExchanged(*scheme, accessor);

  return PtrCouplingScheme(scheme);
}

PtrCouplingScheme CouplingSchemeConfiguration:: createSerialImplicitCouplingScheme
(
  const std::string& accessor ) const
{
  preciceTrace1("createSerialImplicitCouplingScheme()", accessor);
  //assertion1 ( not utils::contained(accessor, _couplingSchemes), accessor );

  m2n::M2N::SharedPointer m2n = _m2nConfig->getM2N (
      _config.participants[0], _config.participants[1] );
  SerialCouplingScheme* scheme = new SerialCouplingScheme (
      _config.maxTime, _config.maxTimesteps, _config.timestepLength,
      _config.validDigits, _config.participants[0], _config.participants[1],
      accessor, m2n, _config.dtMethod, BaseCouplingScheme::Implicit, _config.maxIterations );
  scheme->setCheckPointTimestepInterval(_config.checkpointTimestepInterval);
  scheme->setExtrapolationOrder ( _config.extrapolationOrder );

  addDataToBeExchanged(*scheme, accessor);

  // Add convergence measures
  using boost::get;
  for (size_t i=0; i < _config.convMeasures.size(); i++){
    int dataID = get<0>(_config.convMeasures[i]);
    bool suffices = get<1>(_config.convMeasures[i]);
    std::string neededMesh = get<2>(_config.convMeasures[i]);
    impl::PtrConvergenceMeasure measure = get<3>(_config.convMeasures[i]);
    _meshConfig->addNeededMesh(_config.participants[1],neededMesh);
    checkIfDataIsExchanged(dataID);
    scheme->addConvergenceMeasure(dataID, suffices, measure);
  }

  // Set relaxation parameters
  if (_postProcConfig->getPostProcessing().get() != NULL){
    foreach(std::string& neededMesh, _postProcConfig->getNeededMeshes()){
      _meshConfig->addNeededMesh(_config.participants[1],neededMesh);
    }
    for(const int dataID : _postProcConfig->getPostProcessing()->getDataIDs() ){
      checkIfDataIsExchanged(dataID);
    }
    scheme->setIterationPostProcessing(_postProcConfig->getPostProcessing());
  }
  return PtrCouplingScheme(scheme);
}

PtrCouplingScheme CouplingSchemeConfiguration:: createParallelImplicitCouplingScheme
(
  const std::string& accessor ) const
{
  preciceTrace1("createParallelImplicitCouplingScheme()", accessor);
  assertion1 ( not utils::contained(accessor, _couplingSchemes), accessor );
  m2n::M2N::SharedPointer m2n = _m2nConfig->getM2N(
      _config.participants[0], _config.participants[1] );
  ParallelCouplingScheme* scheme = new ParallelCouplingScheme (
      _config.maxTime, _config.maxTimesteps, _config.timestepLength,
      _config.validDigits, _config.participants[0], _config.participants[1],
      accessor, m2n, _config.dtMethod, BaseCouplingScheme::Implicit, _config.maxIterations );
  scheme->setCheckPointTimestepInterval(_config.checkpointTimestepInterval);
  scheme->setExtrapolationOrder ( _config.extrapolationOrder );

  addDataToBeExchanged(*scheme, accessor);

  // Add convergence measures
  using boost::get;
  for (size_t i=0; i < _config.convMeasures.size(); i++){
    int dataID = get<0>(_config.convMeasures[i]);
    bool suffices = get<1>(_config.convMeasures[i]);
    std::string neededMesh = get<2>(_config.convMeasures[i]);
    impl::PtrConvergenceMeasure measure = get<3>(_config.convMeasures[i]);
    _meshConfig->addNeededMesh(_config.participants[1],neededMesh);
    checkIfDataIsExchanged(dataID);
    scheme->addConvergenceMeasure(dataID, suffices, measure);
  }

  // Set relaxation parameters
  if (_postProcConfig->getPostProcessing().get() != NULL){
    foreach(std::string& neededMesh, _postProcConfig->getNeededMeshes()){
      _meshConfig->addNeededMesh(_config.participants[1],neededMesh);
    }
    for(const int dataID : _postProcConfig->getPostProcessing()->getDataIDs() ){
      checkIfDataIsExchanged(dataID);
    }
    scheme->setIterationPostProcessing(_postProcConfig->getPostProcessing());
  }
  return PtrCouplingScheme(scheme);
}

PtrCouplingScheme CouplingSchemeConfiguration:: createMultiCouplingScheme
(
  const std::string& accessor ) const
{
  preciceTrace1("createMultiCouplingScheme()", accessor);
  assertion1 ( not utils::contained(accessor, _couplingSchemes), accessor );

  BaseCouplingScheme* scheme;

  if(accessor == _config.controller){
    std::vector<m2n::M2N::SharedPointer> m2ns;
    for(const std::string& participant : _config.participants){
      m2ns.push_back(_m2nConfig->getM2N (
          _config.controller, participant ));
    }

    scheme = new MultiCouplingScheme (
        _config.maxTime, _config.maxTimesteps, _config.timestepLength,
        _config.validDigits, accessor, m2ns, _config.dtMethod,
         _config.maxIterations );
    scheme->setCheckPointTimestepInterval(_config.checkpointTimestepInterval);
    scheme->setExtrapolationOrder ( _config.extrapolationOrder );

    MultiCouplingScheme* castedScheme = dynamic_cast<MultiCouplingScheme*>(scheme);
    addMultiDataToBeExchanged(*castedScheme, accessor);
  }
  else{
    m2n::M2N::SharedPointer m2n = _m2nConfig->getM2N (
        accessor, _config.controller );
    scheme = new ParallelCouplingScheme (
        _config.maxTime, _config.maxTimesteps, _config.timestepLength,
        _config.validDigits, accessor, _config.controller,
        accessor, m2n, _config.dtMethod, BaseCouplingScheme::Implicit, _config.maxIterations );
    scheme->setCheckPointTimestepInterval(_config.checkpointTimestepInterval);
    scheme->setExtrapolationOrder ( _config.extrapolationOrder );

    addDataToBeExchanged(*scheme, accessor);
  }

  // Add convergence measures
  using boost::get;
  for (size_t i=0; i < _config.convMeasures.size(); i++){
    int dataID = get<0>(_config.convMeasures[i]);
    bool suffices = get<1>(_config.convMeasures[i]);
    std::string neededMesh = get<2>(_config.convMeasures[i]);
    impl::PtrConvergenceMeasure measure = get<3>(_config.convMeasures[i]);
    _meshConfig->addNeededMesh(_config.controller,neededMesh);
    checkIfDataIsExchanged(dataID);
    scheme->addConvergenceMeasure(dataID, suffices, measure);
  }

  // Set relaxation parameters
  if (_postProcConfig->getPostProcessing().get() != NULL){
    foreach(std::string& neededMesh, _postProcConfig->getNeededMeshes()){
      _meshConfig->addNeededMesh(_config.controller,neededMesh);
    }
    for(const int dataID : _postProcConfig->getPostProcessing()->getDataIDs() ){
      checkIfDataIsExchanged(dataID);
    }

    scheme->setIterationPostProcessing(_postProcConfig->getPostProcessing());
  }
  return PtrCouplingScheme(scheme);
}

constants::TimesteppingMethod
CouplingSchemeConfiguration:: getTimesteppingMethod
(
  const std::string& method ) const
{
  preciceTrace1 ( "getTimesteppingMethod()", method );
  if ( method == VALUE_FIXED ){
    return constants::FIXED_DT;
  }
  else if ( method == VALUE_FIRST_PARTICIPANT ){
    return constants::FIRST_PARTICIPANT_SETS_DT;
  }
//  else if ( method == TagTimestepLength::VALUE_SECOND_PARTICIPANT ){
//    return constants::SECOND_PARTICIPANT_SETS_DT;
//  }
  preciceError ( "getTimesteppingMethod()", "Unknown timestepping method \""
                 << method << "\"!" );
}

void CouplingSchemeConfiguration:: addDataToBeExchanged
(
  BaseCouplingScheme& scheme,
  const std::string&  accessor) const
{
  preciceTrace ( "addDataToBeExchanged()");
  using boost::get;
  for (const Config::Exchange& tuple : _config.exchanges){
    mesh::PtrData data = get<0>(tuple);
    mesh::PtrMesh mesh = get<1>(tuple);
    const std::string& from = get<2>(tuple);
    const std::string& to = get<3>(tuple);

    preciceCheck(to != from,"addDataToBeExchanged()",
        "You cannot define an exchange from and to the same participant");

    if (not(utils::contained(from, _config.participants) || from == _config.controller)){
      throw std::string("Participant \"" + from + "\" is not configured for coupling scheme");
    }

    if (not(utils::contained(to, _config.participants) || to == _config.controller)){
      throw std::string("Participant \"" + to + "\" is not configured for coupling scheme");
    }

    bool initialize = get<4>(tuple);
    if (from == accessor){
      scheme.addDataToSend(data, mesh, initialize);
    }
    else if(to == accessor){
      scheme.addDataToReceive(data, mesh, initialize);
    }
    else{
      assertion(_config.type == VALUE_MULTI);
    }

  }
}

void CouplingSchemeConfiguration:: addMultiDataToBeExchanged
(
  MultiCouplingScheme& scheme,
  const std::string&  accessor) const
{
  preciceTrace ( "addMultiDataToBeExchanged()");
  using boost::get;
  for (const Config::Exchange& tuple : _config.exchanges){
    mesh::PtrData data = get<0>(tuple);
    mesh::PtrMesh mesh = get<1>(tuple);
    const std::string& from = get<2>(tuple);
    const std::string& to  = get<3>(tuple);

    if (not(utils::contained(from, _config.participants) || from == _config.controller)){
      throw std::string("Participant \"" + from + "\" is not configured for coupling scheme");
    }

    if (not(utils::contained(to, _config.participants) || to == _config.controller)){
      throw std::string("Participant \"" + to + "\" is not configured for coupling scheme");
    }

    bool initialize = get<4>(tuple);
    if (from == accessor){
      int index = 0;
      for(const std::string& participant : _config.participants){
        preciceDebug("from: " << from << ", to: " << to << ", participant: " << participant);
        if(to == participant){
          break;
        }
        index++;
      }
      assertion2(index < _config.participants.size(), index, _config.participants.size());
      scheme.addDataToSend(data, mesh, initialize, index);
    }
    else {
      int index = 0;
      for(const std::string& participant : _config.participants){
        preciceDebug("from: " << from << ", to: " << to << ", participant: " << participant);
        if(from == participant){
          break;
        }
        index++;
      }
      assertion2(index < _config.participants.size(), index, _config.participants.size());
      scheme.addDataToReceive(data, mesh, initialize, index);
    }
  }
}

void CouplingSchemeConfiguration:: checkIfDataIsExchanged
(
  int dataID) const
{
  bool hasFound = false;
  using boost::get;
  for (const Config::Exchange& tuple : _config.exchanges){
    mesh::PtrData data = get<0>(tuple);
    if( data->getID()==dataID){
      hasFound = true;
    }
  }
  preciceCheck(hasFound,"checkIfDataIsExchanged()",
        "You need to exchange every data that you use for convergence measures"
        <<" and/or the iteration post-processing");
}







}} // namespace precice, cplscheme
