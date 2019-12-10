#include "MappingConfiguration.hpp"
#include "mapping/NearestNeighborMapping.hpp"
#include "mapping/NearestProjectionMapping.hpp"
#include "mapping/RadialBasisFctMapping.hpp"
#include "mapping/PetRadialBasisFctMapping.hpp"
#include "mapping/impl/BasisFunctions.hpp"
#include "mesh/config/MeshConfiguration.hpp"
#include "xml/XMLTag.hpp"
#include "xml/XMLAttribute.hpp"

namespace precice {
namespace mapping {

MappingConfiguration:: MappingConfiguration
(
  xml::XMLTag &                   parent,
  const mesh::PtrMeshConfiguration& meshConfiguration )
:
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

  _meshConfig(meshConfiguration)
{
  PRECICE_ASSERT(_meshConfig);
  using namespace xml;

 auto attrShapeParam  = XMLAttribute<double>( ATTR_SHAPE_PARAM )
      .setDocumentation("Specific shape parameter for RBF basis function.");
 auto attrSupportRadius  = XMLAttribute<double>( ATTR_SUPPORT_RADIUS )
      .setDocumentation("Support radius of each RBF basis function (global choice).");
  auto attrSolverRtol  = makeXMLAttribute( ATTR_SOLVER_RTOL, 1e-9)
      .setDocumentation("Solver relative tolerance for convergence");
  auto attrXDead = makeXMLAttribute(ATTR_X_DEAD, false)
      .setDocumentation("If set to true, the x axis will be ignored for the mapping");
  auto attrYDead = makeXMLAttribute(ATTR_Y_DEAD, false)
      .setDocumentation("If set to true, the y axis will be ignored for the mapping");
  auto attrZDead = makeXMLAttribute(ATTR_Z_DEAD, false)
      .setDocumentation("If set to true, the z axis will be ignored for the mapping");
  auto attrPolynomial = makeXMLAttribute("polynomial", "separate")
      .setDocumentation("Toggles use of the global polynomial")
      .setOptions({"on", "off", "separate"});

  auto attrPreallocation = makeXMLAttribute("preallocation", "tree")
      .setDocumentation("Sets kind of preallocation for PETSc RBF implementation")
      .setOptions({"estimate", "compute", "off", "save", "tree"});

  XMLTag::Occurrence occ = XMLTag::OCCUR_ARBITRARY;
  std::list<XMLTag> tags;
  {
    XMLTag tag(*this, VALUE_RBF_TPS, occ, TAG);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_RBF_MULTIQUADRICS, occ, TAG);
    tag.addAttribute(attrShapeParam);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_RBF_INV_MULTIQUADRICS, occ, TAG);
    tag.addAttribute(attrShapeParam);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_RBF_VOLUME_SPLINES, occ, TAG);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_RBF_GAUSSIAN, occ, TAG);
    tag.addAttribute(attrShapeParam);
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
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_RBF_CPOLYNOMIAL_C6, occ, TAG);
    tag.addAttribute(attrSupportRadius);
    tags.push_back(tag);
  }
  // ---- Petsc RBF declarations ----
  {
    XMLTag tag(*this, VALUE_PETRBF_TPS, occ, TAG);
    tag.addAttribute(attrSolverRtol);
    tag.addAttribute(attrPolynomial);
    tag.addAttribute(attrPreallocation);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_PETRBF_MULTIQUADRICS, occ, TAG);
    tag.addAttribute(attrShapeParam);
    tag.addAttribute(attrSolverRtol);
    tag.addAttribute(attrPolynomial);
    tag.addAttribute(attrPreallocation);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_PETRBF_INV_MULTIQUADRICS, occ, TAG);
    tag.addAttribute(attrShapeParam);
    tag.addAttribute(attrSolverRtol);
    tag.addAttribute(attrPolynomial);
    tag.addAttribute(attrPreallocation);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_PETRBF_VOLUME_SPLINES, occ, TAG);
    tag.addAttribute(attrSolverRtol);
    tag.addAttribute(attrPolynomial);
    tag.addAttribute(attrPreallocation);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_PETRBF_GAUSSIAN, occ, TAG);
    tag.addAttribute(attrSolverRtol);
    tag.addAttribute(attrShapeParam);
    tag.addAttribute(attrPolynomial);
    tag.addAttribute(attrPreallocation);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_PETRBF_CTPS_C2, occ, TAG);
    tag.addAttribute(attrSolverRtol);
    tag.addAttribute(attrSupportRadius);
    tag.addAttribute(attrPolynomial);
    tag.addAttribute(attrPreallocation);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_PETRBF_CPOLYNOMIAL_C0, occ, TAG);
    tag.addAttribute(attrSolverRtol);
    tag.addAttribute(attrSupportRadius);
    tag.addAttribute(attrPolynomial);
    tag.addAttribute(attrPreallocation);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_PETRBF_CPOLYNOMIAL_C6, occ, TAG);
    tag.addAttribute(attrSolverRtol);
    tag.addAttribute(attrSupportRadius);
    tag.addAttribute(attrPolynomial);
    tag.addAttribute(attrPreallocation);
    tags.push_back(tag);
  }
  // Add tags that only RBF mappings use
  for (XMLTag& tag : tags) {
    tag.addAttribute(attrXDead);
    tag.addAttribute(attrYDead);
    tag.addAttribute(attrZDead);
  }
  {
    XMLTag tag(*this, VALUE_NEAREST_NEIGHBOR, occ, TAG);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_NEAREST_PROJECTION, occ, TAG);
    tags.push_back(tag);
  }
  
  auto attrDirection = XMLAttribute<std::string>( ATTR_DIRECTION)
      .setOptions({ VALUE_WRITE, VALUE_READ });

  XMLAttribute<std::string> attrFromMesh(ATTR_FROM);
  XMLAttribute<std::string> attrToMesh(ATTR_TO);

  auto attrConstraint = XMLAttribute<std::string>(ATTR_CONSTRAINT)
      .setOptions({VALUE_CONSERVATIVE, VALUE_CONSISTENT});

  auto attrTiming = makeXMLAttribute(ATTR_TIMING, VALUE_TIMING_INITIAL)
      .setOptions({VALUE_TIMING_INITIAL, VALUE_TIMING_ON_ADVANCE, VALUE_TIMING_ON_DEMAND});

  // Add tags that all mappings use and add to parent tag
  for (XMLTag & tag : tags) {
    tag.addAttribute(attrDirection);
    tag.addAttribute(attrFromMesh);
    tag.addAttribute(attrToMesh);
    tag.addAttribute(attrConstraint);
    tag.addAttribute(attrTiming);
    parent.addSubtag(tag);
  }
}

void MappingConfiguration:: xmlTagCallback
(
  const xml::ConfigurationContext& context,
  xml::XMLTag& tag )
{
  PRECICE_TRACE(tag.getName());
  if (tag.getNamespace() == TAG){
    std::string dir = tag.getStringAttributeValue(ATTR_DIRECTION);
    std::string fromMesh = tag.getStringAttributeValue(ATTR_FROM);
    std::string toMesh = tag.getStringAttributeValue(ATTR_TO);
    std::string type = tag.getName();
    std::string constraint = tag.getStringAttributeValue(ATTR_CONSTRAINT);
    Timing timing = getTiming(tag.getStringAttributeValue(ATTR_TIMING));
    double shapeParameter = 0.0;
    double supportRadius = 0.0;
    double solverRtol = 1e-9;
    bool xDead = false, yDead = false, zDead = false;
    Polynomial polynomial = Polynomial::ON;
    Preallocation preallocation = Preallocation::TREE;
    
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
    if (tag.hasAttribute("polynomial")) {
      std::string strPolynomial = tag.getStringAttributeValue("polynomial");
      if (strPolynomial == "separate")
        polynomial = Polynomial::SEPARATE;
      else if (strPolynomial == "on")
        polynomial = Polynomial::ON;
      else if (strPolynomial == "off")
        polynomial = Polynomial::OFF;
    }
    if (tag.hasAttribute("preallocation")){
      std::string strPrealloc = tag.getStringAttributeValue("preallocation");
      if (strPrealloc == "estimate")
        preallocation = Preallocation::ESTIMATE;
      else if (strPrealloc == "compute")
        preallocation = Preallocation::COMPUTE;
      else if (strPrealloc == "save")
        preallocation = Preallocation::SAVE;
      else if (strPrealloc == "tree")
        preallocation = Preallocation::TREE;
      else if (strPrealloc == "off")
        preallocation = Preallocation::OFF;
    }     
          
    ConfiguredMapping configuredMapping = createMapping(dir, type, constraint,
                                                        fromMesh, toMesh, timing,
                                                        shapeParameter, supportRadius, solverRtol,
                                                        xDead, yDead, zDead, polynomial, preallocation);
    checkDuplicates ( configuredMapping );
    _mappings.push_back ( configuredMapping );
  }
}

void MappingConfiguration::xmlEndTagCallback(const xml::ConfigurationContext& context, xml::XMLTag& tag)
{
}

const std::vector<MappingConfiguration::ConfiguredMapping>&
MappingConfiguration:: mappings()
{
  return _mappings;
}

MappingConfiguration::ConfiguredMapping MappingConfiguration::createMapping
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
  bool               zDead,
  Polynomial         polynomial,
  Preallocation      preallocation) const
{
  PRECICE_TRACE(direction, type, timing, shapeParameter, supportRadius);
  using namespace mapping;
  ConfiguredMapping configuredMapping;
  mesh::PtrMesh fromMesh(_meshConfig->getMesh(fromMeshName));
  mesh::PtrMesh toMesh(_meshConfig->getMesh(toMeshName));
  PRECICE_CHECK(fromMesh.get() != nullptr, "Mesh \"" << fromMeshName << "\" not defined at creation of mapping!");
  PRECICE_CHECK(toMesh.get() != nullptr, "Mesh \"" << toMeshName << "\" not defined at creation of mapping!");
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
    PRECICE_ERROR("Unknown direction type \"" << direction << "\"!");
  }

  Mapping::Constraint constraintValue;
  if (constraint == VALUE_CONSERVATIVE){
    constraintValue = Mapping::CONSERVATIVE;
  }
  else if (constraint == VALUE_CONSISTENT){
    constraintValue = Mapping::CONSISTENT;
  }
  else {
    PRECICE_ERROR("Unknown mapping constraint \"" << constraint << "\"!");
  }

  # ifndef PRECICE_NO_PETSC
    // for petsc initialization
    int argc = 1;
    char* arg = new char[8];
    strcpy(arg, "precice");
    char** argv = &arg;
  #endif

  configuredMapping.isRBF = true; //will be overwritten for NN and NP just below

  if (type == VALUE_NEAREST_NEIGHBOR){
    configuredMapping.mapping = PtrMapping (
        new NearestNeighborMapping(constraintValue, dimensions) );
    configuredMapping.isRBF = false;
  }
  else if (type == VALUE_NEAREST_PROJECTION){
    configuredMapping.mapping = PtrMapping (
      new NearestProjectionMapping(constraintValue, dimensions) );
    configuredMapping.isRBF = false;
  }
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
    utils::Petsc::initialize(&argc, &argv);
    configuredMapping.mapping = PtrMapping (
      new PetRadialBasisFctMapping<ThinPlateSplines>(constraintValue, dimensions, ThinPlateSplines(),
                                                     xDead, yDead, zDead, solverRtol, polynomial, preallocation) );
  }
  else if (type == VALUE_PETRBF_MULTIQUADRICS){
    utils::Petsc::initialize(&argc, &argv);
    configuredMapping.mapping = PtrMapping (
      new PetRadialBasisFctMapping<Multiquadrics>(constraintValue, dimensions, Multiquadrics(shapeParameter),
                                                  xDead, yDead, zDead, solverRtol, polynomial, preallocation) );
  }
  else if (type == VALUE_PETRBF_INV_MULTIQUADRICS){
    utils::Petsc::initialize(&argc, &argv);
    configuredMapping.mapping = PtrMapping (
      new PetRadialBasisFctMapping<InverseMultiquadrics>(constraintValue, dimensions, InverseMultiquadrics(shapeParameter),
                                                         xDead, yDead, zDead, solverRtol, polynomial, preallocation) );
  }
  else if (type == VALUE_PETRBF_VOLUME_SPLINES){
    utils::Petsc::initialize(&argc, &argv);
    configuredMapping.mapping = PtrMapping (
      new PetRadialBasisFctMapping<VolumeSplines>(constraintValue, dimensions, VolumeSplines(),
                                                  xDead, yDead, zDead, solverRtol, polynomial, preallocation) );
  }
  else if (type == VALUE_PETRBF_GAUSSIAN){
    utils::Petsc::initialize(&argc, &argv);
    configuredMapping.mapping = PtrMapping(
      new PetRadialBasisFctMapping<Gaussian>(constraintValue, dimensions, Gaussian(shapeParameter),
                                             xDead, yDead, zDead, solverRtol, polynomial, preallocation));
  }
  else if (type == VALUE_PETRBF_CTPS_C2){
    utils::Petsc::initialize(&argc, &argv);
    configuredMapping.mapping = PtrMapping (
      new PetRadialBasisFctMapping<CompactThinPlateSplinesC2>(constraintValue, dimensions, CompactThinPlateSplinesC2(supportRadius),
                                                              xDead, yDead, zDead, solverRtol, polynomial, preallocation) );
  }
  else if (type == VALUE_PETRBF_CPOLYNOMIAL_C0){
    utils::Petsc::initialize(&argc, &argv);
    configuredMapping.mapping = PtrMapping (
      new PetRadialBasisFctMapping<CompactPolynomialC0>(constraintValue, dimensions, CompactPolynomialC0(supportRadius),
                                                        xDead, yDead, zDead, solverRtol, polynomial, preallocation) );
  }
  else if (type == VALUE_PETRBF_CPOLYNOMIAL_C6){
    utils::Petsc::initialize(&argc, &argv);
    configuredMapping.mapping = PtrMapping (new PetRadialBasisFctMapping<CompactPolynomialC6>(constraintValue, dimensions, CompactPolynomialC6(supportRadius),
                                                                                              xDead, yDead, zDead, solverRtol, polynomial, preallocation) );
  }
# else
  else if (type.find("petrbf-") == 0) {
    PRECICE_ERROR("PETRBF mappings are not available as preCICE was built without PETSc.");
  }
# endif
  else {
    PRECICE_ERROR("Unknown mapping type!");
  }
  PRECICE_ASSERT(configuredMapping.mapping);
  #ifndef PRECICE_NO_PETSC
    delete[] arg;
  #endif
  return configuredMapping;
}

void MappingConfiguration:: checkDuplicates(const ConfiguredMapping & mapping)
{
  for ( const ConfiguredMapping & configuredMapping : _mappings ) {
    bool sameFromMesh = mapping.fromMesh->getName() == configuredMapping.fromMesh->getName();
    bool sameToMesh = mapping.toMesh->getName() == configuredMapping.toMesh->getName();
    PRECICE_CHECK( !sameFromMesh, "There cannot be two mappings from mesh \""
           << mapping.fromMesh->getName() << "\"" );
    PRECICE_CHECK( !sameToMesh, "There cannot be two mappings to mesh \""
           << mapping.toMesh->getName() << "\"" );
  }
}

MappingConfiguration::Timing MappingConfiguration:: getTiming(const std::string& timing ) const
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
  PRECICE_ERROR("Unknown timing value \"" << timing << "\"!");
}


}} // namespace precice, mapping

