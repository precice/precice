#include "MappingConfiguration.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <list>
#include <memory>
#include <ostream>
#include <string.h>
#include <utility>
#include "logging/LogMacros.hpp"
#include "mapping/Mapping.hpp"
#include "mapping/NearestNeighborMapping.hpp"
#include "mapping/NearestProjectionMapping.hpp"
#include "mapping/PetRadialBasisFctMapping.hpp"
#include "mapping/RadialBasisFctMapping.hpp"
#include "mapping/impl/BasisFunctions.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/config/MeshConfiguration.hpp"
#include "utils/Parallel.hpp"
#include "utils/Petsc.hpp"
#include "utils/assertion.hpp"
#include "xml/ConfigParser.hpp"
#include "xml/XMLAttribute.hpp"
#include "xml/XMLTag.hpp"

namespace precice {
namespace mapping {

MappingConfiguration::MappingConfiguration(
    xml::XMLTag &                     parent,
    const mesh::PtrMeshConfiguration &meshConfiguration)
    : _meshConfig(meshConfiguration)
{
  PRECICE_ASSERT(_meshConfig);
  using namespace xml;

  auto attrShapeParam = XMLAttribute<double>(ATTR_SHAPE_PARAM)
                            .setDocumentation("Specific shape parameter for RBF basis function.");
  auto attrSupportRadius = XMLAttribute<double>(ATTR_SUPPORT_RADIUS)
                               .setDocumentation("Support radius of each RBF basis function (global choice).");
  auto attrSolverRtol = makeXMLAttribute(ATTR_SOLVER_RTOL, 1e-9)
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
  auto attrUseLU = makeXMLAttribute(ATTR_USE_QR, false)
                       .setDocumentation("If set to true, QR decomposition is used to solve the RBF system");

  XMLTag::Occurrence occ = XMLTag::OCCUR_ARBITRARY;
  std::list<XMLTag>  tags;
  {
    XMLTag tag(*this, VALUE_RBF_TPS, occ, TAG);
    tag.setDocumentation("Global radial-basis-function mapping based on the thin plate splines.");
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_RBF_MULTIQUADRICS, occ, TAG);
    tag.setDocumentation("Global radial-basis-function mapping based on the multiquadrics RBF.");
    tag.addAttribute(attrShapeParam);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_RBF_INV_MULTIQUADRICS, occ, TAG);
    tag.setDocumentation("Global radial-basis-function mapping based on the inverse multiquadrics RBF.");
    tag.addAttribute(attrShapeParam);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_RBF_VOLUME_SPLINES, occ, TAG);
    tag.setDocumentation("Global radial-basis-function mapping based on the volume-splines RBF.");
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_RBF_GAUSSIAN, occ, TAG);
    tag.setDocumentation("Local radial-basis-function mapping based on the Gaussian RBF with a cut-off threshold.");
    tag.addAttribute(attrShapeParam);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_RBF_CTPS_C2, occ, TAG);
    tag.setDocumentation("Local radial-basis-function mapping based on the C2-polynomial RBF.");
    tag.addAttribute(attrSupportRadius);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_RBF_CPOLYNOMIAL_C0, occ, TAG);
    tag.setDocumentation("Local radial-basis-function mapping based on the C0-polynomial RBF.");
    tag.addAttribute(attrSupportRadius);
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_RBF_CPOLYNOMIAL_C6, occ, TAG);
    tag.setDocumentation("Local radial-basis-function mapping based on the C6-polynomial RBF.");
    tag.addAttribute(attrSupportRadius);
    tags.push_back(tag);
  }
  // Add tags that only, but all RBF mappings use
  for (XMLTag &tag : tags) {
    tag.addAttribute(attrSolverRtol);
    tag.addAttribute(attrPolynomial);
    tag.addAttribute(attrPreallocation);
    tag.addAttribute(attrXDead);
    tag.addAttribute(attrYDead);
    tag.addAttribute(attrZDead);
    tag.addAttribute(attrUseLU);
  }
  {
    XMLTag tag(*this, VALUE_NEAREST_NEIGHBOR, occ, TAG);
    tag.setDocumentation("Nearest-neighbour mapping which uses a rstar-spacial index tree to index meshes and run nearest-neighbour queries.");
    tags.push_back(tag);
  }
  {
    XMLTag tag(*this, VALUE_NEAREST_PROJECTION, occ, TAG);
    tag.setDocumentation("Nearest-projection mapping which uses a rstar-spacial index tree to index meshes and locate the nearest projections.");
    tags.push_back(tag);
  }

  auto attrDirection = XMLAttribute<std::string>(ATTR_DIRECTION)
                           .setOptions({VALUE_WRITE, VALUE_READ})
                           .setDocumentation("Write mappings map written data prior to communication, thus in the same participant who writes the data. "
                                             "Read mappings map received data after communication, thus in the same participant who reads the data.");

  auto attrFromMesh = XMLAttribute<std::string>(ATTR_FROM)
                          .setDocumentation("The mesh to map the data from.");

  auto attrToMesh = XMLAttribute<std::string>(ATTR_TO)
                        .setDocumentation("The mesh to map the data to.");

  auto attrConstraint = XMLAttribute<std::string>(ATTR_CONSTRAINT)
                            .setDocumentation("Use conservative to conserve the quantity of the data over the interface such as force or mass. Use consistent for normalized quantities such as temperature or pressure.")
                            .setOptions({VALUE_CONSERVATIVE, VALUE_CONSISTENT});

  auto attrTiming = makeXMLAttribute(ATTR_TIMING, VALUE_TIMING_INITIAL)
                        .setDocumentation("This allows to defer the mapping of the data to advance or to a manual call to mapReadDataTo and mapWriteDataFrom.")
                        .setOptions({VALUE_TIMING_INITIAL, VALUE_TIMING_ON_ADVANCE, VALUE_TIMING_ON_DEMAND});

  // Add tags that all mappings use and add to parent tag
  for (XMLTag &tag : tags) {
    tag.addAttribute(attrDirection);
    tag.addAttribute(attrFromMesh);
    tag.addAttribute(attrToMesh);
    tag.addAttribute(attrConstraint);
    tag.addAttribute(attrTiming);
    parent.addSubtag(tag);
  }
}

void MappingConfiguration::xmlTagCallback(
    const xml::ConfigurationContext &context,
    xml::XMLTag &                    tag)
{
  PRECICE_TRACE(tag.getName());
  if (tag.getNamespace() == TAG) {
    std::string   dir            = tag.getStringAttributeValue(ATTR_DIRECTION);
    std::string   fromMesh       = tag.getStringAttributeValue(ATTR_FROM);
    std::string   toMesh         = tag.getStringAttributeValue(ATTR_TO);
    std::string   type           = tag.getName();
    std::string   constraint     = tag.getStringAttributeValue(ATTR_CONSTRAINT);
    Timing        timing         = getTiming(tag.getStringAttributeValue(ATTR_TIMING));
    double        shapeParameter = 0.0;
    double        supportRadius  = 0.0;
    double        solverRtol     = 1e-9;
    bool          xDead = false, yDead = false, zDead = false;
    bool          useLU         = false;
    Polynomial    polynomial    = Polynomial::ON;
    Preallocation preallocation = Preallocation::TREE;

    if (tag.hasAttribute(ATTR_SHAPE_PARAM)) {
      shapeParameter = tag.getDoubleAttributeValue(ATTR_SHAPE_PARAM);
    }
    if (tag.hasAttribute(ATTR_SUPPORT_RADIUS)) {
      supportRadius = tag.getDoubleAttributeValue(ATTR_SUPPORT_RADIUS);
    }
    if (tag.hasAttribute(ATTR_SOLVER_RTOL)) {
      solverRtol = tag.getDoubleAttributeValue(ATTR_SOLVER_RTOL);
    }
    if (tag.hasAttribute(ATTR_X_DEAD)) {
      xDead = tag.getBooleanAttributeValue(ATTR_X_DEAD);
    }
    if (tag.hasAttribute(ATTR_Y_DEAD)) {
      yDead = tag.getBooleanAttributeValue(ATTR_Y_DEAD);
    }
    if (tag.hasAttribute(ATTR_Z_DEAD)) {
      zDead = tag.getBooleanAttributeValue(ATTR_Z_DEAD);
    }
    if (tag.hasAttribute(ATTR_USE_QR)) {
      useLU = tag.getBooleanAttributeValue(ATTR_USE_QR);
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
    if (tag.hasAttribute("preallocation")) {
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

    ConfiguredMapping configuredMapping = createMapping(context,
                                                        dir, type, constraint,
                                                        fromMesh, toMesh, timing,
                                                        shapeParameter, supportRadius, solverRtol,
                                                        xDead, yDead, zDead,
                                                        useLU,
                                                        polynomial, preallocation);
    checkDuplicates(configuredMapping);
    _mappings.push_back(configuredMapping);
  }
}

void MappingConfiguration::xmlEndTagCallback(const xml::ConfigurationContext &context, xml::XMLTag &tag)
{
}

const std::vector<MappingConfiguration::ConfiguredMapping> &
MappingConfiguration::mappings()
{
  return _mappings;
}

MappingConfiguration::ConfiguredMapping MappingConfiguration::createMapping(
    const xml::ConfigurationContext &context,
    const std::string &              direction,
    const std::string &              type,
    const std::string &              constraint,
    const std::string &              fromMeshName,
    const std::string &              toMeshName,
    Timing                           timing,
    double                           shapeParameter,
    double                           supportRadius,
    double                           solverRtol,
    bool                             xDead,
    bool                             yDead,
    bool                             zDead,
    bool                             useLU,
    Polynomial                       polynomial,
    Preallocation                    preallocation) const
{
  PRECICE_TRACE(direction, type, timing, shapeParameter, supportRadius);
  using namespace mapping;
  ConfiguredMapping configuredMapping;
  mesh::PtrMesh     fromMesh(_meshConfig->getMesh(fromMeshName));
  mesh::PtrMesh     toMesh(_meshConfig->getMesh(toMeshName));
  PRECICE_CHECK(fromMesh.get() != nullptr, "Mesh \"" << fromMeshName << "\" was not found while creating a mapping. Please correct the from=\"" << fromMeshName << "\" attribute.");
  PRECICE_CHECK(toMesh.get() != nullptr, "Mesh \"" << toMeshName << "\" was not found while creating a mapping. Please correct the to=\"" << toMeshName << "\" attribute.");
  configuredMapping.fromMesh = fromMesh;
  configuredMapping.toMesh   = toMesh;
  configuredMapping.timing   = timing;
  int dimensions             = fromMesh->getDimensions();

  if (direction == VALUE_WRITE) {
    configuredMapping.direction = WRITE;
  } else if (direction == VALUE_READ) {
    configuredMapping.direction = READ;
  } else {
    PRECICE_ASSERT(false, "Unknown mapping direction type \"" << direction << "\". Please check the documentation for available options.");
  }

  Mapping::Constraint constraintValue;
  if (constraint == VALUE_CONSERVATIVE) {
    constraintValue = Mapping::CONSERVATIVE;
  } else if (constraint == VALUE_CONSISTENT) {
    constraintValue = Mapping::CONSISTENT;
  } else {
    PRECICE_ASSERT(false, "Unknown mapping constraint \"" << constraint << "\". Please check the documentation for available options.");
  }

  if (type == VALUE_NEAREST_NEIGHBOR) {
    configuredMapping.mapping = PtrMapping(
        new NearestNeighborMapping(constraintValue, dimensions));
    configuredMapping.isRBF = false;
    return configuredMapping;
  } else if (type == VALUE_NEAREST_PROJECTION) {
    configuredMapping.mapping = PtrMapping(
        new NearestProjectionMapping(constraintValue, dimensions));
    configuredMapping.isRBF = false;
    return configuredMapping;
  }

  // the mapping is a RBF mapping

  configuredMapping.isRBF = true;
  bool    usePETSc        = false;
  RBFType rbfType         = RBFType::EIGEN;

#ifndef PRECICE_NO_PETSC
  // for petsc initialization
  int   argc = 1;
  char *arg  = new char[8];
  strcpy(arg, "precice");
  char **argv = &arg;
  utils::Petsc::initialize(&argc, &argv, utils::Parallel::current()->comm);
  delete[] arg;
  usePETSc = true;
#endif

  if (usePETSc && (not useLU)) {
    rbfType = RBFType::PETSc;
  } else {
    rbfType = RBFType::EIGEN;
  }

  if (rbfType == RBFType::EIGEN) {
    PRECICE_DEBUG("Eigen RBF is used");
    if (type == VALUE_RBF_TPS) {
      configuredMapping.mapping = PtrMapping(
          new RadialBasisFctMapping<ThinPlateSplines>(constraintValue, dimensions, ThinPlateSplines(), xDead, yDead, zDead));
    } else if (type == VALUE_RBF_MULTIQUADRICS) {
      configuredMapping.mapping = PtrMapping(
          new RadialBasisFctMapping<Multiquadrics>(
              constraintValue, dimensions, Multiquadrics(shapeParameter), xDead, yDead, zDead));
    } else if (type == VALUE_RBF_INV_MULTIQUADRICS) {
      configuredMapping.mapping = PtrMapping(
          new RadialBasisFctMapping<InverseMultiquadrics>(
              constraintValue, dimensions, InverseMultiquadrics(shapeParameter), xDead, yDead, zDead));
    } else if (type == VALUE_RBF_VOLUME_SPLINES) {
      configuredMapping.mapping = PtrMapping(
          new RadialBasisFctMapping<VolumeSplines>(constraintValue, dimensions, VolumeSplines(), xDead, yDead, zDead));
    } else if (type == VALUE_RBF_GAUSSIAN) {
      configuredMapping.mapping = PtrMapping(
          new RadialBasisFctMapping<Gaussian>(
              constraintValue, dimensions, Gaussian(shapeParameter), xDead, yDead, zDead));
    } else if (type == VALUE_RBF_CTPS_C2) {
      configuredMapping.mapping = PtrMapping(
          new RadialBasisFctMapping<CompactThinPlateSplinesC2>(
              constraintValue, dimensions, CompactThinPlateSplinesC2(supportRadius), xDead, yDead, zDead));
    } else if (type == VALUE_RBF_CPOLYNOMIAL_C0) {
      configuredMapping.mapping = PtrMapping(
          new RadialBasisFctMapping<CompactPolynomialC0>(
              constraintValue, dimensions, CompactPolynomialC0(supportRadius), xDead, yDead, zDead));
    } else if (type == VALUE_RBF_CPOLYNOMIAL_C6) {
      configuredMapping.mapping = PtrMapping(
          new RadialBasisFctMapping<CompactPolynomialC6>(
              constraintValue, dimensions, CompactPolynomialC6(supportRadius), xDead, yDead, zDead));
    } else {
      PRECICE_ERROR("Unknown mapping type!");
    }
  }

#ifndef PRECICE_NO_PETSC

  if (rbfType == RBFType::PETSc) {
    PRECICE_DEBUG("PETSc RBF is used.");
    if (type == VALUE_RBF_TPS) {
      configuredMapping.mapping = PtrMapping(
          new PetRadialBasisFctMapping<ThinPlateSplines>(constraintValue, dimensions, ThinPlateSplines(),
                                                         xDead, yDead, zDead, solverRtol, polynomial, preallocation));
    } else if (type == VALUE_RBF_MULTIQUADRICS) {
      configuredMapping.mapping = PtrMapping(
          new PetRadialBasisFctMapping<Multiquadrics>(constraintValue, dimensions, Multiquadrics(shapeParameter),
                                                      xDead, yDead, zDead, solverRtol, polynomial, preallocation));
    } else if (type == VALUE_RBF_INV_MULTIQUADRICS) {
      configuredMapping.mapping = PtrMapping(
          new PetRadialBasisFctMapping<InverseMultiquadrics>(constraintValue, dimensions, InverseMultiquadrics(shapeParameter),
                                                             xDead, yDead, zDead, solverRtol, polynomial, preallocation));
    } else if (type == VALUE_RBF_VOLUME_SPLINES) {
      configuredMapping.mapping = PtrMapping(
          new PetRadialBasisFctMapping<VolumeSplines>(constraintValue, dimensions, VolumeSplines(),
                                                      xDead, yDead, zDead, solverRtol, polynomial, preallocation));
    } else if (type == VALUE_RBF_GAUSSIAN) {
      configuredMapping.mapping = PtrMapping(
          new PetRadialBasisFctMapping<Gaussian>(constraintValue, dimensions, Gaussian(shapeParameter),
                                                 xDead, yDead, zDead, solverRtol, polynomial, preallocation));
    } else if (type == VALUE_RBF_CTPS_C2) {
      configuredMapping.mapping = PtrMapping(
          new PetRadialBasisFctMapping<CompactThinPlateSplinesC2>(constraintValue, dimensions, CompactThinPlateSplinesC2(supportRadius),
                                                                  xDead, yDead, zDead, solverRtol, polynomial, preallocation));
    } else if (type == VALUE_RBF_CPOLYNOMIAL_C0) {
      configuredMapping.mapping = PtrMapping(
          new PetRadialBasisFctMapping<CompactPolynomialC0>(constraintValue, dimensions, CompactPolynomialC0(supportRadius),
                                                            xDead, yDead, zDead, solverRtol, polynomial, preallocation));
    } else if (type == VALUE_RBF_CPOLYNOMIAL_C6) {
      configuredMapping.mapping = PtrMapping(
          new PetRadialBasisFctMapping<CompactPolynomialC6>(constraintValue, dimensions, CompactPolynomialC6(supportRadius),
                                                            xDead, yDead, zDead, solverRtol, polynomial, preallocation));
    } else {
      PRECICE_ERROR("Unknown mapping type!");
    }
  }

#endif

  PRECICE_ASSERT(configuredMapping.mapping);
  return configuredMapping;
}

void MappingConfiguration::checkDuplicates(const ConfiguredMapping &mapping)
{
  for (const ConfiguredMapping &configuredMapping : _mappings) {
    bool sameToMesh   = mapping.toMesh->getName() == configuredMapping.toMesh->getName();
    bool sameFromMesh = mapping.fromMesh->getName() == configuredMapping.fromMesh->getName();
    bool sameMapping  = sameToMesh && sameFromMesh;
    PRECICE_CHECK(!sameMapping, "There cannot be two mappings from mesh \""
                                    << mapping.fromMesh->getName() << "\" to mesh \""
                                    << mapping.toMesh->getName() << "\". "
                                    << "Please remove one of the duplicated meshes. ");
  }
}

MappingConfiguration::Timing MappingConfiguration::getTiming(const std::string &timing) const
{
  if (timing == VALUE_TIMING_INITIAL) {
    return INITIAL;
  } else if (timing == VALUE_TIMING_ON_ADVANCE) {
    return ON_ADVANCE;
  } else if (timing == VALUE_TIMING_ON_DEMAND) {
    return ON_DEMAND;
  }
  PRECICE_ASSERT(false, "Unknown timing value \"" << timing << "\". Please check the documentation for available options.");
}

} // namespace mapping
} // namespace precice
