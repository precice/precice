// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "CouplingSchemeConfiguration.hpp"
#include "cplscheme/config/PostProcessingConfiguration.hpp"
#include "cplscheme/ExplicitCouplingScheme.hpp"
#include "cplscheme/ImplicitCouplingScheme.hpp"
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
#include "com/Communication.hpp"
#include "com/config/CommunicationConfiguration.hpp"
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
  const com::PtrCommunicationConfiguration& comConfig )
:
  TAG("coupling-scheme"),
  TAG_PARTICIPANTS("participants"),
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
  ATTR_SUFFICES("suffices"),
  VALUE_EXPLICIT("explicit"),
  VALUE_IMPLICIT("implicit"),
  VALUE_UNCOUPLED("uncoupled"),
  VALUE_FIXED("fixed"),
  VALUE_FIRST_PARTICIPANT("first-participant"),
  _config(),
  _meshConfig(meshConfig),
  _comConfig(comConfig),
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
    XMLTag tag(*this, VALUE_EXPLICIT, occ, TAG);
    doc = "Explicit coupling scheme according to conventional serial";
    doc += " staggered prcedure (CSS).";
    tag.setDocumentation(doc);
    addTypespecifcSubtags(VALUE_EXPLICIT, tag);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_IMPLICIT, occ, TAG);
    doc = "Implicit coupling scheme according to block Gauss-Seidel iterations.";
    doc += " Improved implicit iterations are achieved by using a post-processing.";
    tag.setDocumentation(doc);
    addTypespecifcSubtags(VALUE_IMPLICIT, tag);
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
    _config.participant = tag.getStringAttributeValue(ATTR_FIRST);
    _config.secondParticipant = tag.getStringAttributeValue(ATTR_SECOND);
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
    assertion(_config.type == VALUE_IMPLICIT);
    addAbsoluteConvergenceMeasure(dataName, meshName, limit, suffices);
  }
  else if ( tag.getName() == TAG_REL_CONV_MEASURE ) {
    std::string dataName = tag.getStringAttributeValue (ATTR_DATA);
    std::string meshName = tag.getStringAttributeValue(ATTR_MESH);
    double limit = tag.getDoubleAttributeValue(ATTR_LIMIT);
    bool suffices = tag.getBooleanAttributeValue(ATTR_SUFFICES);
    assertion(_config.type == VALUE_IMPLICIT);
    addRelativeConvergenceMeasure(dataName, meshName, limit, suffices);
  }
  else if ( tag.getName() == TAG_RES_REL_CONV_MEASURE ) {
    std::string dataName = tag.getStringAttributeValue (ATTR_DATA);
    std::string meshName = tag.getStringAttributeValue(ATTR_MESH);
    double limit = tag.getDoubleAttributeValue(ATTR_LIMIT);
    bool suffices = tag.getBooleanAttributeValue(ATTR_SUFFICES);
    assertion(_config.type == VALUE_IMPLICIT);
    addResidualRelativeConvergenceMeasure(dataName, meshName, limit, suffices);
  }
  else if ( tag.getName() == TAG_MIN_ITER_CONV_MEASURE ) {
    std::string dataName = tag.getStringAttributeValue(ATTR_DATA);
    std::string meshName = tag.getStringAttributeValue(ATTR_MESH);
    int minIterations = tag.getIntAttributeValue(ATTR_MIN_ITERATIONS);
    bool suffices = tag.getBooleanAttributeValue(ATTR_SUFFICES);
    assertion(_config.type == VALUE_IMPLICIT);
    addMinIterationConvergenceMeasure(dataName, meshName, minIterations, suffices);
  }
  else if (tag.getName() == TAG_EXCHANGE){
    assertion(_config.type != VALUE_UNCOUPLED);
    std::string nameData = tag.getStringAttributeValue(ATTR_DATA);
    std::string nameMesh = tag.getStringAttributeValue(ATTR_MESH);
    std::string nameParticipant = tag.getStringAttributeValue(ATTR_FROM);
    bool initialize = tag.getBooleanAttributeValue(ATTR_INITIALIZE);
    mesh::PtrData exchangeData;
    foreach (mesh::PtrMesh mesh, _meshConfig->meshes()){
      if ( mesh->getName() == nameMesh ) {
        foreach ( mesh::PtrData data, mesh->data() ) {
          if ( data->getName() == nameData ) {
            exchangeData = data;
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
    _config.exchanges.push_back(boost::make_tuple(exchangeData,
                                nameParticipant, initialize));
  }
  else if (tag.getName() == TAG_MAX_ITERATIONS){
    assertion(_config.type == VALUE_IMPLICIT);
    _config.maxIterations = tag.getIntAttributeValue(ATTR_VALUE);
  }
  else if (tag.getName() == TAG_EXTRAPOLATION){
    assertion(_config.type == VALUE_IMPLICIT);
    _config.extrapolationOrder = tag.getIntAttributeValue(ATTR_VALUE);
  }
}

void CouplingSchemeConfiguration:: xmlEndTagCallback
(
  utils::XMLTag& tag )
{
  preciceTrace1("xmlEndTagCallback()", tag.getFullName());
  if (tag.getNamespace() == TAG){
    if (_config.type == VALUE_EXPLICIT){
      std::string accessor(_config.participant);
      PtrCouplingScheme scheme = createExplicitCouplingScheme(accessor);
      addCouplingScheme(scheme, accessor);
      //_couplingSchemes[accessor] = scheme;
      accessor = _config.secondParticipant;
      scheme = createExplicitCouplingScheme(accessor);
      addCouplingScheme(scheme, accessor);
      //_couplingSchemes[accessor] = scheme;
      _config = Config();
    }
    else if (_config.type == VALUE_IMPLICIT){
      std::string accessor(_config.participant);
      PtrCouplingScheme scheme = createImplicitCouplingScheme(accessor);
      addCouplingScheme(scheme, accessor);
      //_couplingSchemes[accessor] = scheme;
      accessor = _config.secondParticipant;
      scheme = createImplicitCouplingScheme(accessor);
      addCouplingScheme(scheme, accessor);
      //_couplingSchemes[accessor] = scheme;
      _config = Config();
    }
    else if (_config.type == VALUE_UNCOUPLED){
      assertion(false);
    }
    else {
      assertion(false);
    }
  }
}

void CouplingSchemeConfiguration:: addCouplingScheme
(
  PtrCouplingScheme  cplScheme,
  const std::string& participantName )
{
  preciceTrace1 ( "addCouplingScheme()", participantName );
  if (_couplingSchemes.find(participantName) != _couplingSchemes.end()){
    preciceDebug("Coupling scheme exists already for participant");
    if (_couplingSchemeCompositions.find(participantName)
        != _couplingSchemeCompositions.end())
    {
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
  addTransientLimitTags(tag);
  _config.type = type;
  //_config.name = name;

  addTagCheckpoint(tag);

  if (type == VALUE_EXPLICIT){
    addTagParticipants(tag);
    addTagExchange(tag);
  }
  else if ( type == VALUE_IMPLICIT ) {
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
  XMLAttribute<std::string> participant(ATTR_FROM);
  tagExchange.addAttribute(participant);
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
  _postProcConfig = PtrPostProcessingConfiguration(
                    new PostProcessingConfiguration(tag, _meshConfig));
}

void CouplingSchemeConfiguration:: addAbsoluteConvergenceMeasure
(
  const std::string& dataName,
  const std::string& meshName,
  double             limit,
  bool               suffices )
{
  impl::PtrConvergenceMeasure measure(new impl::AbsoluteConvergenceMeasure(limit));
  int dataID = getData(dataName, meshName)->getID();
  _config.convMeasures.push_back(boost::make_tuple(dataID, suffices, measure));
}

void CouplingSchemeConfiguration:: addRelativeConvergenceMeasure
(
  const std::string& dataName,
  const std::string& meshName,
  double             limit,
  bool               suffices )
{
  impl::PtrConvergenceMeasure measure (
      new impl::RelativeConvergenceMeasure(limit) );
  int dataID = getData(dataName, meshName)->getID();
  _config.convMeasures.push_back(boost::make_tuple(dataID, suffices, measure));
}

void CouplingSchemeConfiguration:: addResidualRelativeConvergenceMeasure
(
  const std::string& dataName,
  const std::string& meshName,
  double             limit,
  bool               suffices )
{
  impl::PtrConvergenceMeasure measure (
      new impl::ResidualRelativeConvergenceMeasure(limit) );
  int dataID = getData(dataName, meshName)->getID();
  _config.convMeasures.push_back(boost::make_tuple(dataID, suffices, measure));
}

void CouplingSchemeConfiguration:: addMinIterationConvergenceMeasure
(
  const std::string& dataName,
  const std::string& meshName,
  int                minIterations,
  bool               suffices )
{
  impl::PtrConvergenceMeasure measure (
      new impl::MinIterationConvergenceMeasure(minIterations) );
  int dataID = getData(dataName, meshName)->getID();
  _config.convMeasures.push_back(boost::make_tuple(dataID, suffices, measure));

}

mesh::PtrData CouplingSchemeConfiguration:: getData
(
  const std::string& dataName,
  const std::string& meshName ) const
{
  foreach ( mesh::PtrMesh mesh, _meshConfig->meshes() ){
    if ( meshName == mesh->getName() ){
      foreach ( mesh::PtrData data, mesh->data() ){
        if ( dataName == data->getName() ){
          return data;
        }
      }
    }
  }
  preciceError ( "getData()", "Data \"" << dataName << "\" used by mesh \""
                 << meshName << "\" is not configured!" );
}


PtrCouplingScheme CouplingSchemeConfiguration:: createExplicitCouplingScheme
(
  const std::string& accessor ) const
{
  //assertion ( not utils::contained(accessor, _couplingSchemes) );
  com::PtrCommunication com = _comConfig->getCommunication (
      _config.participant, _config.secondParticipant );
  ExplicitCouplingScheme* scheme = new ExplicitCouplingScheme (
      _config.maxTime, _config.maxTimesteps, _config.timestepLength,
      _config.validDigits, _config.participant, _config.secondParticipant,
      accessor, com, _config.dtMethod );
  scheme->setCheckointTimestepInterval ( _config.checkpointTimestepInterval );

  addDataToBeExchanged(*scheme, accessor);

  return PtrCouplingScheme(scheme);
}

PtrCouplingScheme CouplingSchemeConfiguration:: createImplicitCouplingScheme
(
  const std::string& accessor ) const
{
  //assertion1 ( not utils::contained(accessor, _couplingSchemes), accessor );
  com::PtrCommunication com = _comConfig->getCommunication (
      _config.participant, _config.secondParticipant );
  ImplicitCouplingScheme* scheme = new ImplicitCouplingScheme (
      _config.maxTime, _config.maxTimesteps, _config.timestepLength,
      _config.validDigits, _config.participant, _config.secondParticipant,
      accessor, com, _config.maxIterations, _config.dtMethod );
  scheme->setCheckointTimestepInterval(_config.checkpointTimestepInterval);
  scheme->setExtrapolationOrder ( _config.extrapolationOrder );

  addDataToBeExchanged(*scheme, accessor);
  // Add data to be sent and received
//  using boost::get;
//  foreach ( const Config::Exchange & tuple, _config.exchanges ) {
//    mesh::PtrData data = get<0> ( tuple );
//    const std::string & from = get<1> ( tuple );
//    bool initialize = get<2> ( tuple );
//    if ( from == accessor ) {
//      scheme->addDataToSend ( data, initialize );
//    }
//    else {
//      scheme->addDataToReceive ( data, initialize );
//    }
//  }

  // Add convergence measures
  using boost::get;
  for (size_t i=0; i < _config.convMeasures.size(); i++){
    int dataID = get<0>(_config.convMeasures[i]);
    bool suffices = get<1>(_config.convMeasures[i]);
    impl::PtrConvergenceMeasure measure = get<2>(_config.convMeasures[i]);
    scheme->addConvergenceMeasure(dataID, suffices, measure);
  }

  // Set relaxation parameters
  if (_postProcConfig->getPostProcessing().get() != NULL){
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
  using boost::get;
  foreach (const Config::Exchange& tuple, _config.exchanges){
    mesh::PtrData data = get<0>(tuple);
    const std::string& from = get<1>(tuple);
    if (from.compare(_config.participant) && from.compare(_config.secondParticipant)){
      throw std::string("Participant \"" + from + "\" is not configured for coupling scheme");
    }

    bool initialize = get<2>(tuple);
    if (from == accessor){
      scheme.addDataToSend(data, initialize);
    }
    else {
      scheme.addDataToReceive(data, initialize);
    }
  }
}

}} // namespace precice, cplscheme
