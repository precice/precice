// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "ActionConfiguration.hpp"
#include "utils/xml/XMLAttribute.hpp"
#include "utils/xml/ValidatorEquals.hpp"
#include "utils/xml/ValidatorOr.hpp"
#include "precice/impl/ModifyCoordinatesDataAction.hpp"
#include "precice/impl/ScaleByAreaAction.hpp"
#include "precice/impl/ScaleByDtAction.hpp"
#include "precice/impl/ComputeCurvatureDataAction.hpp"
#include "precice/impl/BalanceVertexPositionAction.hpp"

namespace precice {
namespace config {

tarch::logging::Log ActionConfiguration::
  _log ( "precice::config::ActionConfiguration" );

const std::string & ActionConfiguration:: getTag ()
{
  static std::string tag = "action";
  return tag;
}

ActionConfiguration:: ActionConfiguration
(
  utils::ptr_vector<impl::MeshContext> & usedMeshes )
:
  TAG_SOURCE_DATA ( "source-data" ),
  TAG_TARGET_DATA ( "target-data" ),
  TAG_MESH ( "mesh" ),
  TAG_CONVERGENCE_TOLERANCE ( "convergence-tolerance" ),
  TAG_MAX_ITERATIONS ( "max-iterations" ),
  ATTR_TYPE ( "type" ),
  ATTR_TIMING ( "timing" ),
  ATTR_NAME ( "name" ),
  ATTR_VALUE ( "value" ),
  VALUE_REGULAR_PRIOR ( "regular-prior" ),
  VALUE_REGULAR_POST ( "regular-post" ),
  VALUE_ON_EXCHANGE_PRIOR ( "on-exchange-prior" ),
  VALUE_ON_EXCHANGE_POST ( "on-exchange-post" ),
  VALUE_DIVIDE_BY_AREA ( "divide-by-area" ),
  VALUE_MULTIPLY_BY_AREA ( "multiply-by-area" ),
  VALUE_SCALE_BY_COMPUTED_DT ( "scale-by-computed-dt" ),
  VALUE_SCALE_BY_COMPUTED_DT_PART ( "scale-by-computed-dt-part" ),
  //VALUE_SET_AS_COORDINATES ( "set-as-coordinates" ),
  VALUE_ADD_TO_COORDINATES ( "add-to-coordinates" ),
  VALUE_SUBTRACT_FROM_COORDINATES ( "subtract-from-coordinates" ),
  VALUE_COMPUTE_CURVATURE ( "compute-curvature" ),
  VALUE_BALANCE_VERTEX_POSITIONS ( "balance-vertex-positions" ),
  _isValid ( false ),
  _usedMeshes ( usedMeshes ),
  _action ()
//  _xmlTag ( getTag(), utils::XMLTag::OCCUR_ONCE )
{}

bool ActionConfiguration:: parseSubtag
(
  tarch::irr::io::IrrXMLReader * xmlReader )
{
  preciceTrace ( "parseSubtag()" );
  using namespace utils;
  _configured = Configured();

  XMLTag tag ( getTag(), utils::XMLTag::OCCUR_ONCE );

  XMLAttribute<std::string> attrType ( ATTR_TYPE );
  ValidatorEquals<std::string> validDivideByArea ( VALUE_DIVIDE_BY_AREA );
  ValidatorEquals<std::string> validMultiplyByArea ( VALUE_MULTIPLY_BY_AREA );
  ValidatorEquals<std::string> validScaleByComputedDt ( VALUE_SCALE_BY_COMPUTED_DT );
  ValidatorEquals<std::string> validScaleByComputedDtPart ( VALUE_SCALE_BY_COMPUTED_DT_PART );
//  ValidatorEquals<std::string> validSetAsCoordinates ( VALUE_SET_AS_COORDINATES );
  ValidatorEquals<std::string> validAddToCoordinates ( VALUE_ADD_TO_COORDINATES );
  ValidatorEquals<std::string> validSubtractFromCoordinates ( VALUE_SUBTRACT_FROM_COORDINATES );
  ValidatorEquals<std::string> validComputeCurvature ( VALUE_COMPUTE_CURVATURE );
  ValidatorEquals<std::string> validBalanceVertexPos ( VALUE_BALANCE_VERTEX_POSITIONS );
  attrType.setValidator ( validDivideByArea || validMultiplyByArea || validScaleByComputedDt
      || validScaleByComputedDtPart /*|| validSetAsCoordinates*/ || validAddToCoordinates
      || validSubtractFromCoordinates || validComputeCurvature || validBalanceVertexPos );
  tag.addAttribute ( attrType );

  XMLAttribute<std::string> attrTiming ( ATTR_TIMING );
  ValidatorEquals<std::string> validRegularPrior ( VALUE_REGULAR_PRIOR );
  ValidatorEquals<std::string> validRegularPost ( VALUE_REGULAR_POST );
  ValidatorEquals<std::string> validOnExchangePrior ( VALUE_ON_EXCHANGE_PRIOR );
  ValidatorEquals<std::string> validOnExchangePost ( VALUE_ON_EXCHANGE_POST );
  attrTiming.setValidator ( validRegularPrior || validRegularPost ||
                            validOnExchangePrior || validOnExchangePost );
  tag.addAttribute ( attrTiming );

  XMLTag tagMesh ( TAG_MESH, XMLTag::OCCUR_ONCE );
  XMLAttribute<std::string> attrName ( ATTR_NAME );
  tagMesh.addAttribute ( attrName );
  tag.addSubtag ( tagMesh );

  _isValid = tag.parse ( xmlReader, *this );
  if ( _isValid ){
    createAction ();
  }
  return _isValid;
}

bool ActionConfiguration:: xmlTagCallback
(
  utils::XMLTag &                callingTag,
  tarch::irr::io::IrrXMLReader * xmlReader )
{
  preciceTrace1 ( "xmlTagCallback()", callingTag.getName() );
  if ( callingTag.getName() == getTag() ){
    _configured = Configured();
    _configured.type = callingTag.getStringAttribute(ATTR_TYPE).getValue();
    _configured.timing = callingTag.getStringAttribute(ATTR_TIMING).getValue();
    addSubtags ( callingTag, _configured.type );
  }
  else if ( callingTag.getName() == TAG_SOURCE_DATA ){
    _configured.sourceData = callingTag.getStringAttribute(ATTR_NAME).getValue();
  }
  else if ( callingTag.getName() == TAG_TARGET_DATA ){
    _configured.targetData = callingTag.getStringAttribute(ATTR_NAME).getValue();
  }
  else if ( callingTag.getName() == TAG_MESH ){
    _configured.mesh = callingTag.getStringAttribute(ATTR_NAME).getValue();
  }
  else if ( callingTag.getName() == TAG_CONVERGENCE_TOLERANCE ){
    _configured.convergenceTolerance =
        callingTag.getDoubleAttribute(ATTR_VALUE).getValue();
  }
  else if ( callingTag.getName() == TAG_MAX_ITERATIONS ){
    _configured.maxIterations = callingTag.getIntAttribute(ATTR_VALUE).getValue();
  }
  return true;
}

bool ActionConfiguration:: xmlEndTagCallback
(
  utils::XMLTag &                callingTag,
  tarch::irr::io::IrrXMLReader * xmlReader )
{
  return true;
}

bool ActionConfiguration:: isValid () const
{
  return _isValid;
}

const impl::PtrAbstractDataAction & ActionConfiguration:: getAction () const
{
  assertion ( _isValid );
  assertion ( _action.use_count() > 0 );
  return _action;
}

void ActionConfiguration:: addSubtags
(
  utils::XMLTag &     callingTag,
  const std::string & type )
{
  preciceTrace1 ( "addSubtags()", callingTag.getName() );
  assertion ( type != std::string("") );
  using utils::XMLTag;
  using utils::XMLAttribute;
  XMLTag tagSourceData ( TAG_SOURCE_DATA, XMLTag::OCCUR_ONCE );
  XMLTag tagTargetData ( TAG_TARGET_DATA, XMLTag::OCCUR_ONCE );

  XMLAttribute<std::string> attrName ( ATTR_NAME );
  tagSourceData.addAttribute ( attrName );
  tagTargetData.addAttribute ( attrName );

  if ( _configured.type == VALUE_ADD_TO_COORDINATES ){
    callingTag.addSubtag ( tagSourceData );
  }
//  else if ( callingTag.getName() == VALUE_SET_AS_COORDINATES ){
//    callingTag.addSubtag ( tagSourceData );
//  }
  else if ( _configured.type == VALUE_SUBTRACT_FROM_COORDINATES ){
    callingTag.addSubtag ( tagSourceData );
  }
  else if ( _configured.type == VALUE_MULTIPLY_BY_AREA ){
    callingTag.addSubtag ( tagTargetData );
  }
  else if ( _configured.type == VALUE_DIVIDE_BY_AREA ){
    callingTag.addSubtag ( tagTargetData );
  }
  else if ( _configured.type == VALUE_SCALE_BY_COMPUTED_DT ){
    callingTag.addSubtag ( tagSourceData );
    callingTag.addSubtag ( tagTargetData );
  }
  else if ( _configured.type == VALUE_SCALE_BY_COMPUTED_DT_PART ){
    callingTag.addSubtag ( tagSourceData );
    callingTag.addSubtag ( tagTargetData );
  }
  else if ( _configured.type == VALUE_COMPUTE_CURVATURE ){
    callingTag.addSubtag ( tagTargetData );
  }
  else if ( _configured.type == VALUE_BALANCE_VERTEX_POSITIONS ){
    XMLTag tagConvergenceTol ( TAG_CONVERGENCE_TOLERANCE, XMLTag::OCCUR_ONCE );
    XMLTag tagMaxIterations ( TAG_MAX_ITERATIONS, XMLTag::OCCUR_ONCE );
    XMLAttribute<double> attrDoubleValue ( ATTR_VALUE );
    XMLAttribute<int> attrIntValue ( ATTR_VALUE );
    tagConvergenceTol.addAttribute ( attrDoubleValue );
    tagMaxIterations.addAttribute ( attrIntValue );
    callingTag.addSubtag ( tagConvergenceTol );
    callingTag.addSubtag ( tagMaxIterations );
  }
  else {
    preciceError ( "addSubtags()", "Unknown action type \""
                   << _configured.type << "\"!" );
  }
}

void ActionConfiguration:: createAction ()
{
  preciceTrace ( "createAction()" );

  assertion ( _configured.type != std::string("") );
  impl::AbstractDataAction::TimingConstants timing = getTiming();

  // Determine data and mesh
  int sourceDataID = -1;
  int targetDataID = -1;
  impl::MeshContext * meshContext = NULL;
  foreach ( impl::MeshContext & context, _usedMeshes ){
    if ( context.mesh->getName() == _configured.mesh ){
      meshContext = & context;
      foreach ( const mesh::PtrData & data, context.mesh->data() ){
        if ( data->getName() == _configured.sourceData ){
          sourceDataID = data->getID ();
        }
        if ( data->getName() == _configured.targetData ){
          targetDataID = data->getID ();
        }
      }
    }
  }
  preciceCheck ( meshContext != NULL, "xmlEndTagCallback()",
                 "Data action uses mesh \"" << _configured.mesh
                 << "\" which is not used!" );

  if ( _configured.type == VALUE_ADD_TO_COORDINATES ){
    typedef impl::ModifyCoordinatesDataAction Action;
    Action::ModeConstants mode = Action::ADD_TO_COORDINATES_MODE;
    _action = impl::PtrAbstractDataAction (
        new Action(timing, sourceDataID, *meshContext, mode) );
  }
//  else if ( callingTag.getName() == VALUE_SET_AS_COORDINATES ){
//
//  }
  else if ( _configured.type == VALUE_SUBTRACT_FROM_COORDINATES ){
    typedef impl::ModifyCoordinatesDataAction Action;
    Action::ModeConstants mode = Action::SUBTRACT_FROM_COORDINATES_MODE;
    _action = impl::PtrAbstractDataAction (
        new Action(timing, sourceDataID, *meshContext, mode) );
  }
  else if ( _configured.type == VALUE_MULTIPLY_BY_AREA ){
    _action = impl::PtrAbstractDataAction (
        new impl::ScaleByAreaAction ( timing, targetDataID,
        meshContext->mesh, impl::ScaleByAreaAction::SCALING_MULTIPLY_BY_AREA) );
  }
  else if ( _configured.type == VALUE_DIVIDE_BY_AREA ){
    _action = impl::PtrAbstractDataAction (
        new impl::ScaleByAreaAction ( timing, targetDataID,
        meshContext->mesh, impl::ScaleByAreaAction::SCALING_DIVIDE_BY_AREA) );
  }
  else if ( _configured.type == VALUE_SCALE_BY_COMPUTED_DT ){
    _action = impl::PtrAbstractDataAction (
        new impl::ScaleByDtAction ( timing, sourceDataID, targetDataID,
        meshContext->mesh, impl::ScaleByDtAction::SCALING_BY_COMPUTED_TIMESTEP) );
  }
  else if ( _configured.type == VALUE_SCALE_BY_COMPUTED_DT_PART ){
    _action = impl::PtrAbstractDataAction (
        new impl::ScaleByDtAction ( timing, sourceDataID, targetDataID,
        meshContext->mesh, impl::ScaleByDtAction::SCALING_BY_COMPUTED_TIMESTEP_PART) );
  }
  else if ( _configured.type == VALUE_COMPUTE_CURVATURE ){
    _action = impl::PtrAbstractDataAction (
        new impl::ComputeCurvatureDataAction ( timing, sourceDataID,
        meshContext->mesh ) );
  }
  else if ( _configured.type == VALUE_BALANCE_VERTEX_POSITIONS ){
    _action = impl::PtrAbstractDataAction (
        new impl::BalanceVertexPositionAction(timing, meshContext->mesh,
        _configured.convergenceTolerance, _configured.maxIterations) );
  }
  else {
    preciceError ( "xmlEndtagCallback()", "Unknown action type \""
                   << _configured.type << "\"!" );
  }
}

impl::AbstractDataAction::TimingConstants ActionConfiguration:: getTiming () const
{
  preciceTrace1 ( "getTiming()", _configured.timing );
  impl::AbstractDataAction::TimingConstants timing;
  if ( _configured.timing == VALUE_REGULAR_PRIOR ){
    timing = impl::AbstractDataAction::ALWAYS_PRIOR;
  }
  else if ( _configured.timing == VALUE_REGULAR_POST ){
    timing = impl::AbstractDataAction::ALWAYS_POST;
  }
  else if ( _configured.timing == VALUE_ON_EXCHANGE_PRIOR ){
    timing = impl::AbstractDataAction::ON_EXCHANGE_PRIOR;
  }
  else if ( _configured.timing == VALUE_ON_EXCHANGE_POST ){
    timing = impl::AbstractDataAction::ON_EXCHANGE_POST;
  }
  else {
    preciceError ( "getTiming()", "Unknown action timing \""
                   <<  _configured.timing << "\"!" );
  }
  return timing;
}

}} // namespace precice, config

