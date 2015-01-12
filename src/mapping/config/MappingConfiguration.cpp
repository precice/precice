// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "MappingConfiguration.hpp"
#include "mapping/NearestNeighborMapping.hpp"
#include "mapping/NearestProjectionMapping.hpp"
#include "mapping/RadialBasisFctMapping.hpp"
#include "mapping/PetRadialBasisFctMapping.hpp"
#include "mesh/config/MeshConfiguration.hpp"
#include "utils/Globals.hpp"
#include "utils/xml/XMLTag.hpp"
#include "utils/xml/XMLAttribute.hpp"
#include "utils/xml/ValidatorEquals.hpp"
#include "utils/xml/ValidatorOr.hpp"

namespace precice {
namespace mapping {

tarch::logging::Log MappingConfiguration::
    _log ( "precice::config::MappingConfiguration" );

//const std::string& MappingConfiguration:: getTag ()
//{
//  static std::string tag("mapping");
//  return tag;
//}

MappingConfiguration:: MappingConfiguration
(
  utils::XMLTag&                    parent,
  const mesh::PtrMeshConfiguration& meshConfiguration )
:
  TAG("mapping"),
  ATTR_DIRECTION("direction"),
  ATTR_FROM("from"),
  ATTR_TO("to"),
  ATTR_TIMING("timing"),
  ATTR_TYPE("type"),
  ATTR_CONSTRAINT("constraint"),
  ATTR_SHAPE_PARAM("shape-parameter"),
  ATTR_SUPPORT_RADIUS("support-radius"),
  ATTR_SOLVER_RTOL("solver-rtol"),
  ATTR_X_DEAD("x-dead"),
  ATTR_Y_DEAD("y-dead"),
  ATTR_Z_DEAD("z-dead"),
  VALUE_WRITE("write"),
  VALUE_READ("read"),
  VALUE_CONSISTENT("consistent"),
  VALUE_CONSERVATIVE("conservative"),
  VALUE_NEAREST_NEIGHBOR("nearest-neighbor"),
  VALUE_NEAREST_PROJECTION("nearest-projection"),
  VALUE_RBF_TPS("rbf-thin-plate-splines"),
  VALUE_RBF_MULTIQUADRICS("rbf-multiquadrics"),
  VALUE_RBF_INV_MULTIQUADRICS("rbf-inverse-multiquadrics"),
  VALUE_RBF_VOLUME_SPLINES("rbf-volume-splines"),
  VALUE_RBF_GAUSSIAN("rbf-gaussian"),
  VALUE_RBF_CTPS_C2("rbf-compact-tps-c2"),
  VALUE_RBF_CPOLYNOMIAL_C0("rbf-compact-polynomial-c0"),
  VALUE_RBF_CPOLYNOMIAL_C6("rbf-compact-polynomial-c6"),

  VALUE_PETRBF_TPS("petrbf-thin-plate-splines"),
  VALUE_PETRBF_MULTIQUADRICS("petrbf-multiquadrics"),
  VALUE_PETRBF_INV_MULTIQUADRICS("petrbf-inverse-multiquadrics"),
  VALUE_PETRBF_VOLUME_SPLINES("petrbf-volume-splines"),
  VALUE_PETRBF_GAUSSIAN("petrbf-gaussian"),
  VALUE_PETRBF_CTPS_C2("petrbf-compact-tps-c2"),
  VALUE_PETRBF_CPOLYNOMIAL_C0("petrbf-compact-polynomial-c0"),
  VALUE_PETRBF_CPOLYNOMIAL_C6("petrbf-compact-polynomial-c6"),

  VALUE_TIMING_INITIAL("initial"),
  VALUE_TIMING_ON_ADVANCE("onadvance"),
  VALUE_TIMING_ON_DEMAND("ondemand"),
  _meshConfig(meshConfiguration),
  //_isValid(false),
  _mappings()
{
  assertion (_meshConfig.use_count() > 0);
  using namespace utils;

  XMLAttribute<double> attrShapeParam ( ATTR_SHAPE_PARAM );
  XMLAttribute<double> attrSupportRadius ( ATTR_SUPPORT_RADIUS );
  XMLAttribute<double> attrSolverRtol ( ATTR_SOLVER_RTOL );
  XMLAttribute<bool> attrXDead(ATTR_X_DEAD);
  std::string doc = "If set to true, the x axis will be ignored for the mapping";
  attrXDead.setDocumentation(doc);
  attrXDead.setDefaultValue(false);
  XMLAttribute<bool> attrYDead(ATTR_Y_DEAD);
  doc = "If set to true, the y axis will be ignored for the mapping";
  attrYDead.setDocumentation(doc);
  attrYDead.setDefaultValue(false);
  XMLAttribute<bool> attrZDead(ATTR_Z_DEAD);
  doc = "If set to true, the z axis will be ignored for the mapping";
  attrZDead.setDocumentation(doc);
  attrZDead.setDefaultValue(false);



  XMLTag::Occurrence occ = XMLTag::OCCUR_ARBITRARY;
  std::list<XMLTag> tags;
  {
    XMLTag tag(*this, VALUE_NEAREST_NEIGHBOR, occ, TAG);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_NEAREST_PROJECTION, occ, TAG);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_RBF_TPS, occ, TAG);
    tag.addAttribute(attrXDead);
    tag.addAttribute(attrYDead);
    tag.addAttribute(attrZDead);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_RBF_MULTIQUADRICS, occ, TAG);
    tag.addAttribute(attrShapeParam);
    tag.addAttribute(attrXDead);
    tag.addAttribute(attrYDead);
    tag.addAttribute(attrZDead);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_RBF_INV_MULTIQUADRICS, occ, TAG);
    tag.addAttribute(attrShapeParam);
    tag.addAttribute(attrXDead);
    tag.addAttribute(attrYDead);
    tag.addAttribute(attrZDead);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_RBF_VOLUME_SPLINES, occ, TAG);
    tag.addAttribute(attrXDead);
    tag.addAttribute(attrYDead);
    tag.addAttribute(attrZDead);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_RBF_GAUSSIAN, occ, TAG);
    tag.addAttribute(attrShapeParam);
    tag.addAttribute(attrXDead);
    tag.addAttribute(attrYDead);
    tag.addAttribute(attrZDead);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_RBF_CTPS_C2, occ, TAG);
    tag.addAttribute(attrSupportRadius);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_RBF_CPOLYNOMIAL_C0, occ, TAG);
    tag.addAttribute(attrSupportRadius);
    tag.addAttribute(attrXDead);
    tag.addAttribute(attrYDead);
    tag.addAttribute(attrZDead);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_RBF_CPOLYNOMIAL_C6, occ, TAG);
    tag.addAttribute(attrSupportRadius);
    tag.addAttribute(attrXDead);
    tag.addAttribute(attrYDead);
    tag.addAttribute(attrZDead);
    tags.push_back(tag);
  }
  // ---- Petsc RBF declarations ----
  {
    XMLTag tag(*this, VALUE_PETRBF_TPS, occ, TAG);
    tag.addAttribute(attrSolverRtol);
    tag.addAttribute(attrXDead);
    tag.addAttribute(attrYDead);
    tag.addAttribute(attrZDead);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_PETRBF_MULTIQUADRICS, occ, TAG);
    tag.addAttribute(attrShapeParam);
    tag.addAttribute(attrSolverRtol);
    tag.addAttribute(attrXDead);
    tag.addAttribute(attrYDead);
    tag.addAttribute(attrZDead);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_PETRBF_INV_MULTIQUADRICS, occ, TAG);
    tag.addAttribute(attrShapeParam);
    tag.addAttribute(attrSolverRtol);
    tag.addAttribute(attrXDead);
    tag.addAttribute(attrYDead);
    tag.addAttribute(attrZDead);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_PETRBF_VOLUME_SPLINES, occ, TAG);
    tag.addAttribute(attrSolverRtol);
    tag.addAttribute(attrXDead);
    tag.addAttribute(attrYDead);
    tag.addAttribute(attrZDead);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_PETRBF_GAUSSIAN, occ, TAG);
    tag.addAttribute(attrSolverRtol);
    tag.addAttribute(attrShapeParam);
    tag.addAttribute(attrXDead);
    tag.addAttribute(attrYDead);
    tag.addAttribute(attrZDead);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_PETRBF_CTPS_C2, occ, TAG);
    tag.addAttribute(attrSolverRtol);
    tag.addAttribute(attrSupportRadius);
    tag.addAttribute(attrXDead);
    tag.addAttribute(attrYDead);
    tag.addAttribute(attrZDead);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_PETRBF_CPOLYNOMIAL_C0, occ, TAG);
    tag.addAttribute(attrSolverRtol);
    tag.addAttribute(attrSupportRadius);
    tag.addAttribute(attrXDead);
    tag.addAttribute(attrYDead);
    tag.addAttribute(attrZDead);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_PETRBF_CPOLYNOMIAL_C6, occ, TAG);
    tag.addAttribute(attrSolverRtol);
    tag.addAttribute(attrSupportRadius);
    tag.addAttribute(attrXDead);
    tag.addAttribute(attrYDead);
    tag.addAttribute(attrZDead);
    tags.push_back(tag);
  }

  XMLAttribute<std::string> attrDirection ( ATTR_DIRECTION );
  ValidatorEquals<std::string> validDirectionWrite ( VALUE_WRITE );
  ValidatorEquals<std::string> validDirectionRead ( VALUE_READ );
  attrDirection.setValidator ( validDirectionWrite || validDirectionRead );


  XMLAttribute<std::string> attrFromMesh(ATTR_FROM);
  XMLAttribute<std::string> attrToMesh(ATTR_TO);

//  XMLAttribute<std::string> attrType ( ATTR_TYPE );
  typedef ValidatorEquals<std::string> ValidString;
//  ValidString validNN ( VALUE_NEAREST_NEIGHBOR );
//  ValidString validNP ( VALUE_NEAREST_PROJECTION );
//  ValidString validRBFTPS ( VALUE_RBF_TPS );
//  ValidString validRBFMultiquadrics ( VALUE_RBF_MULTIQUADRICS );
//  ValidString validRBFInvMultiquadrics ( VALUE_RBF_INV_MULTIQUADRICS );
//  ValidString validRBFVolumeSplines ( VALUE_RBF_VOLUME_SPLINES );
//  ValidString validRBFGaussian ( VALUE_RBF_GAUSSIAN );
//  ValidString validRBFCTPSC2 ( VALUE_RBF_CTPS_C2 );
//  ValidString validRBFCPolynomialC0 ( VALUE_RBF_CPOLYNOMIAL_C0 );
//  ValidString validRBFCPolynomialC6 ( VALUE_RBF_CPOLYNOMIAL_C6 );
//  attrType.setValidator ( validNN || validNP || validRBFTPS || validRBFMultiquadrics
//      || validRBFInvMultiquadrics || validRBFVolumeSplines || validRBFGaussian
//      || validRBFCTPSC2 || validRBFCPolynomialC0 || validRBFCPolynomialC6 );
//  tagMapping.addAttribute ( attrType );

  XMLAttribute<std::string> attrConstraint(ATTR_CONSTRAINT);
  ValidString validConservative(VALUE_CONSERVATIVE);
  ValidString validConsistent(VALUE_CONSISTENT);
  attrConstraint.setValidator(validConservative || validConsistent);

  XMLAttribute<std::string> attrTiming(ATTR_TIMING);
  attrTiming.setDefaultValue(VALUE_TIMING_INITIAL);
  ValidString validInitial(VALUE_TIMING_INITIAL);
  ValidString validOnAdvance(VALUE_TIMING_ON_ADVANCE);
  ValidString validOnDemand(VALUE_TIMING_ON_DEMAND);
  attrTiming.setValidator(validInitial || validOnAdvance || validOnDemand);

  foreach (XMLTag& tag, tags){
    tag.addAttribute(attrDirection);
    tag.addAttribute(attrFromMesh);
    tag.addAttribute(attrToMesh);
    tag.addAttribute(attrConstraint);
    tag.addAttribute(attrTiming);
    //tag.addAttribute(attrIncremental);
    parent.addSubtag(tag);
  }
}


//bool MappingConfiguration:: parseSubtag
//(
//  utils::XMLTag::XMLReader* xmlReader )
//{
//  using utils::XMLTag;
//  using utils::XMLAttribute;
//  using utils::ValidatorEquals;
//
//  XMLTag tagMapping ( TAG, XMLTag::OCCUR_ONCE );
//
//  XMLAttribute<std::string> attrDirection ( ATTR_DIRECTION );
//  ValidatorEquals<std::string> validDirectionWrite ( VALUE_WRITE );
//  ValidatorEquals<std::string> validDirectionRead ( VALUE_READ );
//  attrDirection.setValidator ( validDirectionWrite || validDirectionRead );
//  tagMapping.addAttribute ( attrDirection );
//
//  XMLAttribute<std::string> attrMesh ( ATTR_MESH );
//  attrMesh.setDefaultValue ( "" );
//  tagMapping.addAttribute ( attrMesh );
//
//  XMLAttribute<std::string> attrType ( ATTR_TYPE );
//  typedef ValidatorEquals<std::string> ValidString;
//  ValidString validNN ( VALUE_NEAREST_NEIGHBOR );
//  ValidString validNP ( VALUE_NEAREST_PROJECTION );
//  ValidString validRBFTPS ( VALUE_RBF_TPS );
//  ValidString validRBFMultiquadrics ( VALUE_RBF_MULTIQUADRICS );
//  ValidString validRBFInvMultiquadrics ( VALUE_RBF_INV_MULTIQUADRICS );
//  ValidString validRBFVolumeSplines ( VALUE_RBF_VOLUME_SPLINES );
//  ValidString validRBFGaussian ( VALUE_RBF_GAUSSIAN );
//  ValidString validRBFCTPSC2 ( VALUE_RBF_CTPS_C2 );
//  ValidString validRBFCPolynomialC0 ( VALUE_RBF_CPOLYNOMIAL_C0 );
//  ValidString validRBFCPolynomialC6 ( VALUE_RBF_CPOLYNOMIAL_C6 );
//  attrType.setValidator ( validNN || validNP || validRBFTPS || validRBFMultiquadrics
//      || validRBFInvMultiquadrics || validRBFVolumeSplines || validRBFGaussian
//      || validRBFCTPSC2 || validRBFCPolynomialC0 || validRBFCPolynomialC6 );
//  tagMapping.addAttribute ( attrType );
//
//  XMLAttribute<std::string> attrConstraint(ATTR_CONSTRAINT);
//  ValidString validConservative(VALUE_CONSERVATIVE);
//  ValidString validConsistent(VALUE_CONSISTENT);
//  attrConstraint.setValidator(validConservative || validConsistent);
//  tagMapping.addAttribute(attrConstraint);
//
//  XMLAttribute<bool> attrStationary ( ATTR_STATIONARY );
//  attrStationary.setDefaultValue ( false );
//  tagMapping.addAttribute ( attrStationary );
//
//  XMLAttribute<double> attrShapeParam ( ATTR_SHAPE_PARAM );
//  attrShapeParam.setDefaultValue ( 0.0 );
//  tagMapping.addAttribute ( attrShapeParam );
//
//  XMLAttribute<double> attrSupportRadius ( ATTR_SUPPORT_RADIUS );
//  attrSupportRadius.setDefaultValue ( std::numeric_limits<double>::max() );
//  tagMapping.addAttribute ( attrSupportRadius );
//
//  _isValid = tagMapping.parse ( xmlReader, *this );
//  return _isValid;
//}

void MappingConfiguration:: xmlTagCallback
(
  utils::XMLTag& tag )
{
  preciceTrace1("xmlTagCallback()", tag.getName());
  if (tag.getNamespace() == TAG){
    std::string dir = tag.getStringAttributeValue(ATTR_DIRECTION);
    std::string fromMesh = tag.getStringAttributeValue(ATTR_FROM);
    std::string toMesh = tag.getStringAttributeValue(ATTR_TO);
    std::string type = tag.getName(); //StringAttributeValue(ATTR_TYPE);
    std::string constraint = tag.getStringAttributeValue(ATTR_CONSTRAINT);
    Timing timing = getTiming(tag.getStringAttributeValue(ATTR_TIMING));
    double shapeParameter = 0.0;
    double supportRadius = 0.0;
    double solverRtol = 1e-9;
    bool xDead = false;
    bool yDead = false;
    bool zDead = false;
    if (tag.hasAttribute(ATTR_SHAPE_PARAM)){
      shapeParameter = tag.getDoubleAttributeValue(ATTR_SHAPE_PARAM);
    }
    if (tag.hasAttribute(ATTR_SUPPORT_RADIUS)){
      supportRadius = tag.getDoubleAttributeValue(ATTR_SUPPORT_RADIUS);
    }
    if (tag.hasAttribute(ATTR_SOLVER_RTOL)){
      solverRtol = tag.getDoubleAttributeValue(ATTR_SOLVER_RTOL);
    }
    if (tag.hasAttribute(ATTR_X_DEAD)){
      xDead = tag.getBooleanAttributeValue(ATTR_X_DEAD);
    }
    if (tag.hasAttribute(ATTR_Y_DEAD)){
      yDead = tag.getBooleanAttributeValue(ATTR_Y_DEAD);
    }
    if (tag.hasAttribute(ATTR_Z_DEAD)){
      zDead = tag.getBooleanAttributeValue(ATTR_Z_DEAD);
    }
        
    ConfiguredMapping configuredMapping = createMapping(dir, type, constraint,
      fromMesh, toMesh, timing, shapeParameter, supportRadius, solverRtol,
      xDead, yDead, zDead);
    checkDuplicates ( configuredMapping );
    _mappings.push_back ( configuredMapping );
  }
}

void MappingConfiguration:: xmlEndTagCallback
(
  utils::XMLTag& tag )
{
}


//bool MappingConfiguration:: isValid () const
//{
//  return _isValid;
//}

const std::vector<MappingConfiguration::ConfiguredMapping>&
MappingConfiguration:: mappings()
{
  return _mappings;
}

void MappingConfiguration:: addMapping
(
  const PtrMapping&    mapping,
  const mesh::PtrMesh& fromMesh,
  const mesh::PtrMesh& toMesh,
  Direction            direction,
  Timing               timing )
{
  preciceTrace3("addMapping()", fromMesh, direction, timing);
  //assertion ( !((timing == INITIALLY) && isIncremental) );
  ConfiguredMapping configuredMapping;
  configuredMapping.mapping = mapping;
  configuredMapping.fromMesh = fromMesh;
  configuredMapping.toMesh = toMesh;
  configuredMapping.direction = direction;
  //configuredMapping.isIncremental = isIncremental;
  configuredMapping.timing = timing;
  checkDuplicates ( configuredMapping );
  _mappings.push_back ( configuredMapping );
}

MappingConfiguration::ConfiguredMapping MappingConfiguration:: createMapping
(
  const std::string& direction,
  const std::string& type,
  const std::string& constraint,
  const std::string& fromMeshName,
  const std::string& toMeshName,
  Timing             timing,
  double             shapeParameter,
  double             supportRadius,
  double             solverRtol,
  bool               xDead,
  bool               yDead,
  bool               zDead) const
{
  preciceTrace5("createMapping()", direction, type, timing,
                shapeParameter, supportRadius);
  using namespace mapping;
  ConfiguredMapping configuredMapping;
  mesh::PtrMesh fromMesh(_meshConfig->getMesh(fromMeshName));
  mesh::PtrMesh toMesh(_meshConfig->getMesh(toMeshName));
  preciceCheck(fromMesh.get() != NULL, "createMapping()",
               "Mesh \"" << fromMeshName << "\" not defined at creation of mapping!");
  preciceCheck(toMesh.get() != NULL, "createMapping()",
               "Mesh \"" << toMeshName << "\" not defined at creation of mapping!");
  configuredMapping.fromMesh = fromMesh;
  configuredMapping.toMesh = toMesh;
  configuredMapping.timing = timing;
  int dimensions = fromMesh->getDimensions();

  if (direction == VALUE_WRITE){
    configuredMapping.direction =  WRITE;
  }
  else if (direction == VALUE_READ){
    configuredMapping.direction = READ;
  }
  else {
    preciceError("createMapping()", "Unknown direction type \"" << direction << "\"!");
  }

  Mapping::Constraint constraintValue;
  if (constraint == VALUE_CONSERVATIVE){
    constraintValue = Mapping::CONSERVATIVE;
  }
  else if (constraint == VALUE_CONSISTENT){
    constraintValue = Mapping::CONSISTENT;
  }
  else {
    preciceError("createMapping()",
                 "Unknown mapping constraint \"" << constraint << "\"!");
  }

  if (type == VALUE_NEAREST_NEIGHBOR){
    configuredMapping.mapping = PtrMapping (
        new NearestNeighborMapping(constraintValue, dimensions) );
  }
  else if (type == VALUE_NEAREST_PROJECTION){
//    preciceCheck ( direction == VALUE_WRITE, "createMapping()",
//                   "Attribute \"direction\" has to have value \"" << VALUE_WRITE
//                   << "\" for a " << "mapping of type \""
//                   << VALUE_CONSERVATIVE_NEAREST_PROJECTION << "\"!" );
    configuredMapping.mapping = PtrMapping (
        new NearestProjectionMapping(constraintValue, dimensions) );
  }
//  else if (type == VALUE_CONSISTENT_NEAREST_PROJECTION){
//    preciceCheck ( direction == VALUE_READ, "createMapping()",
//                   "Attribute \"direction\" has to have value \"" << VALUE_READ
//                   << "\" for a " << "mapping of type \""
//                   << VALUE_CONSISTENT_NEAREST_PROJECTION << "\"!" );
//    configuredMapping.mapping = PtrMapping (
//        new MappingConsistentNearestProjection() );
//  }
  else if (type == VALUE_RBF_TPS){
    configuredMapping.mapping = PtrMapping (
      new RadialBasisFctMapping<ThinPlateSplines>(constraintValue, dimensions, ThinPlateSplines(),
            xDead, yDead, zDead));
  }
  else if (type == VALUE_RBF_MULTIQUADRICS){
    configuredMapping.mapping = PtrMapping (
      new RadialBasisFctMapping<Multiquadrics>(
        constraintValue, dimensions, Multiquadrics(shapeParameter),
        xDead, yDead, zDead ));
  }
  else if (type == VALUE_RBF_INV_MULTIQUADRICS){
    configuredMapping.mapping = PtrMapping (
      new RadialBasisFctMapping<InverseMultiquadrics>(
        constraintValue, dimensions, InverseMultiquadrics(shapeParameter),
        xDead, yDead, zDead ));
  }
  else if (type == VALUE_RBF_VOLUME_SPLINES){
    configuredMapping.mapping = PtrMapping (
      new RadialBasisFctMapping<VolumeSplines>(constraintValue, dimensions, VolumeSplines(),
      xDead, yDead, zDead ));
  }
  else if (type == VALUE_RBF_GAUSSIAN){
    configuredMapping.mapping = PtrMapping(
        new RadialBasisFctMapping<Gaussian>(
          constraintValue, dimensions, Gaussian(shapeParameter),
          xDead, yDead, zDead));
  }
  else if (type == VALUE_RBF_CTPS_C2){
    configuredMapping.mapping = PtrMapping (
      new RadialBasisFctMapping<CompactThinPlateSplinesC2>(
        constraintValue, dimensions, CompactThinPlateSplinesC2(supportRadius),
        xDead, yDead, zDead ));
  }
  else if (type == VALUE_RBF_CPOLYNOMIAL_C0){
    configuredMapping.mapping = PtrMapping (
      new RadialBasisFctMapping<CompactPolynomialC0>(
        constraintValue, dimensions, CompactPolynomialC0(supportRadius),
        xDead, yDead, zDead ));
  }
  else if (type == VALUE_RBF_CPOLYNOMIAL_C6){
    configuredMapping.mapping = PtrMapping (
      new RadialBasisFctMapping<CompactPolynomialC6>(
        constraintValue, dimensions, CompactPolynomialC6(supportRadius),
        xDead, yDead, zDead ));
  }
# ifndef PRECICE_NO_PETSC
  else if (type == VALUE_PETRBF_TPS){
    configuredMapping.mapping = PtrMapping (
      new PetRadialBasisFctMapping<ThinPlateSplines>(constraintValue, dimensions, ThinPlateSplines(),
                                                     xDead, yDead, zDead, solverRtol) );
  }
  else if (type == VALUE_PETRBF_MULTIQUADRICS){
    configuredMapping.mapping = PtrMapping (
      new PetRadialBasisFctMapping<Multiquadrics>(
        constraintValue, dimensions, Multiquadrics(shapeParameter),
        xDead, yDead, zDead, solverRtol) );
  }
  else if (type == VALUE_PETRBF_INV_MULTIQUADRICS){
    configuredMapping.mapping = PtrMapping (
      new PetRadialBasisFctMapping<InverseMultiquadrics>(
        constraintValue, dimensions, InverseMultiquadrics(shapeParameter),
        xDead, yDead, zDead, solverRtol) );
  }
  else if (type == VALUE_PETRBF_VOLUME_SPLINES){
    configuredMapping.mapping = PtrMapping (
      new PetRadialBasisFctMapping<VolumeSplines>(constraintValue, dimensions, VolumeSplines(),
                                                  xDead, yDead, zDead, solverRtol) );
  }
  else if (type == VALUE_PETRBF_GAUSSIAN){
    configuredMapping.mapping = PtrMapping(
        new PetRadialBasisFctMapping<Gaussian>(
          constraintValue, dimensions, Gaussian(shapeParameter),
          xDead, yDead, zDead, solverRtol));
  }
  else if (type == VALUE_PETRBF_CTPS_C2){
    configuredMapping.mapping = PtrMapping (
      new PetRadialBasisFctMapping<CompactThinPlateSplinesC2>(
        constraintValue, dimensions, CompactThinPlateSplinesC2(supportRadius),
        xDead, yDead, zDead, solverRtol) );
  }
  else if (type == VALUE_PETRBF_CPOLYNOMIAL_C0){
    configuredMapping.mapping = PtrMapping (
      new PetRadialBasisFctMapping<CompactPolynomialC0>(
        constraintValue, dimensions, CompactPolynomialC0(supportRadius),
        xDead, yDead, zDead, solverRtol) );
  }
  else if (type == VALUE_PETRBF_CPOLYNOMIAL_C6){
    configuredMapping.mapping = PtrMapping (
      new PetRadialBasisFctMapping<CompactPolynomialC6>(
        constraintValue, dimensions, CompactPolynomialC6(supportRadius),
        xDead, yDead, zDead, solverRtol) );
  }
# endif
  else {
    preciceError ( "getMapping()", "Unknown mapping type!" );
  }
  assertion ( configuredMapping.mapping.use_count() > 0 );
  return configuredMapping;
}

void MappingConfiguration:: checkDuplicates
(
  const ConfiguredMapping & mapping )
{
  foreach ( const ConfiguredMapping & configuredMapping, _mappings ) {
    bool sameFromMesh = mapping.fromMesh->getName() == configuredMapping.fromMesh->getName();
    bool sameToMesh = mapping.toMesh->getName() == configuredMapping.toMesh->getName();
    preciceCheck ( !sameFromMesh, "checkDuplicates()",
                   "There cannot be two mappings from mesh \""
                   << mapping.fromMesh->getName() << "\"" );
    preciceCheck ( !sameToMesh, "checkDuplicates()",
                       "There cannot be two mappings to mesh \""
                       << mapping.toMesh->getName() << "\"" );
  }
}

MappingConfiguration::Timing MappingConfiguration:: getTiming(
  const std::string& timing ) const
{
  if (timing == VALUE_TIMING_INITIAL){
    return INITIAL;
  }
  else if (timing == VALUE_TIMING_ON_ADVANCE){
    return ON_ADVANCE;
  }
  else if (timing == VALUE_TIMING_ON_DEMAND){
    return ON_DEMAND;
  }
  preciceError("getTiming()", "Unknown timing value \"" << timing << "\"!");
}


}} // namespace precice, mapping
