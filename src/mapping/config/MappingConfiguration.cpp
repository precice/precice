#include "MappingConfiguration.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <cstring>
#include <list>
#include <memory>
#include <ostream>
#include <utility>
#include <variant>
#include "logging/LogMacros.hpp"
#include "mapping/LinearCellInterpolationMapping.hpp"
#include "mapping/Mapping.hpp"
#include "mapping/NearestNeighborGradientMapping.hpp"
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
namespace precice::mapping {

namespace {

template <typename Listener, typename TagStorage>
void createTag(Listener &              listener,
               const std::string &     name,
               xml::XMLTag::Occurrence occurrence,
               const std::string &     xmlNamespace,
               TagStorage &            storage,
               std::string             documentation)
{
  xml::XMLTag tag(listener, name, occurrence, xmlNamespace);
  tag.setDocumentation(documentation);
  storage.push_back(tag);
}

void addSubtagsToParents(std::list<xml::XMLTag> &subtags,
                         std::list<xml::XMLTag> &parents)
{
  for (auto &p : parents) {
    std::for_each(subtags.begin(), subtags.end(), [&p](auto &s) { p.addSubtag(s); });
  }
}

using variant_t = std::variant<xml::XMLAttribute<double>, xml::XMLAttribute<std::string>, xml::XMLAttribute<bool>>;
template <typename TagStorage>
void addAttributes(TagStorage &storage, const std::vector<variant_t> &attributes)
{
  for (auto &s : storage) {
    for (auto &a : attributes)
      std::visit([&s](auto &&arg) { s.addAttribute(arg); }, a);
  }
}
} // namespace

MappingConfiguration::MappingConfiguration(
    xml::XMLTag &              parent,
    mesh::PtrMeshConfiguration meshConfiguration)
    : _meshConfig(std::move(meshConfiguration))
{
  PRECICE_ASSERT(_meshConfig);
  using namespace xml;

  // First, we create the available tags
  XMLTag::Occurrence occ = XMLTag::OCCUR_ARBITRARY;
  std::list<XMLTag>  projectionTags, rbfDirectTags, rbfIterativeTags;
  createTag(*this, TYPE_NEAREST_NEIGHBOR, occ, TAG, projectionTags, "Nearest-neighbour mapping which uses a rstar-spacial index tree to index meshes and run nearest-neighbour queries.");
  createTag(*this, TYPE_NEAREST_PROJECTION, occ, TAG, projectionTags, "Nearest-projection mapping which uses a rstar-spacial index tree to index meshes and locate the nearest projections.");
  createTag(*this, TYPE_NEAREST_NEIGHBOR_GRADIENT, occ, TAG, projectionTags, "Nearest-neighbor-gradient mapping which uses nearest-neighbor mapping with an additional linear approximation using gradient data.");
  createTag(*this, TYPE_LINEAR_CELL_INTERPOLATION, occ, TAG, projectionTags, "Linear cell interpolation mapping which uses a rstar-spacial index tree to index meshes and locate the nearest cell. Only supports 2D meshes.");
  createTag(*this, TYPE_RBF_GLOBAL_DIRECT, occ, TAG, rbfDirectTags, "Radial-basis-function mapping using a direct solver with a gather-scatter parallelism.");
  createTag(*this, TYPE_RBF_GLOBAL_ITERATIVE, occ, TAG, rbfIterativeTags, "Radial-basis-function mapping using an iterative solver with a distributed parallelism.");

  // List of all attributes with corresponding documentation
  auto attrDirection = XMLAttribute<std::string>(ATTR_DIRECTION)
                           .setOptions({DIRECTION_WRITE, DIRECTION_READ})
                           .setDocumentation("Write mappings map written data prior to communication, thus in the same participant who writes the data. "
                                             "Read mappings map received data after communication, thus in the same participant who reads the data.");

  auto attrFromMesh = XMLAttribute<std::string>(ATTR_FROM)
                          .setDocumentation("The mesh to map the data from.");

  auto attrToMesh = XMLAttribute<std::string>(ATTR_TO)
                        .setDocumentation("The mesh to map the data to.");

  auto attrConstraint = XMLAttribute<std::string>(ATTR_CONSTRAINT)
                            .setDocumentation("Use conservative to conserve the nodal sum of the data over the interface (needed e.g. for force mapping).  Use consistent for normalized quantities such as temperature or pressure. Use scaled-consistent-surface or scaled-consistent-volume for normalized quantities where conservation of integral values (surface or volume) is needed (e.g. velocities when the mass flow rate needs to be conserved). Mesh connectivity is required to use scaled-consistent.")
                            .setOptions({CONSTRAINT_CONSERVATIVE, CONSTRAINT_CONSISTENT, CONSTRAINT_SCALED_CONSISTENT_SURFACE, CONSTRAINT_SCALED_CONSISTENT_VOLUME});
  auto attrXDead = makeXMLAttribute(ATTR_X_DEAD, false)
                       .setDocumentation("If set to true, the x axis will be ignored for the mapping");
  auto attrYDead = makeXMLAttribute(ATTR_Y_DEAD, false)
                       .setDocumentation("If set to true, the y axis will be ignored for the mapping");
  auto attrZDead = makeXMLAttribute(ATTR_Z_DEAD, false)
                       .setDocumentation("If set to true, the z axis will be ignored for the mapping");
  auto attrPolynomial = makeXMLAttribute(ATTR_POLYNOMIAL, POLYNOMIAL_SEPARATE)
                            .setDocumentation("Toggles use of the global polynomial")
                            .setOptions({POLYNOMIAL_ON, POLYNOMIAL_OFF, POLYNOMIAL_SEPARATE});

  auto attrSolverRtol = makeXMLAttribute(ATTR_SOLVER_RTOL, 1e-9)
                            .setDocumentation("Solver relative tolerance for convergence");
  auto attrPreallocation = makeXMLAttribute(ATTR_PREALLOCATION, PREALLOCATION_TREE)
                               .setDocumentation("Sets kind of preallocation for PETSc RBF implementation")
                               .setOptions({PREALLOCATION_ESTIMATE, PREALLOCATION_COMPUTE, PREALLOCATION_OFF, PREALLOCATION_SAVE, PREALLOCATION_TREE});

  // Add the relevant attributes to the relevant tags
  addAttributes(projectionTags, {attrFromMesh, attrToMesh, attrDirection, attrConstraint});
  addAttributes(rbfDirectTags, {attrFromMesh, attrToMesh, attrDirection, attrConstraint, attrPolynomial, attrXDead, attrYDead, attrZDead});
  addAttributes(rbfIterativeTags, {attrFromMesh, attrToMesh, attrDirection, attrConstraint, attrPolynomial, attrXDead, attrYDead, attrZDead, attrSolverRtol, attrPreallocation});

  // Now we take care of the subtag basis function
  // First, we have the tags using a support radius
  XMLTag::Occurrence once = XMLTag::OCCUR_NOT_OR_ONCE;
  std::list<XMLTag>  supportRadiusRBF;
  createTag(*this, RBF_CPOLYNOMIAL_C0, once, SUBTAG_BASIS_FUNCTION, supportRadiusRBF, "Wendland C0 function");
  createTag(*this, RBF_CPOLYNOMIAL_C2, once, SUBTAG_BASIS_FUNCTION, supportRadiusRBF, "Wendland C2 function");
  createTag(*this, RBF_CPOLYNOMIAL_C4, once, SUBTAG_BASIS_FUNCTION, supportRadiusRBF, "Wendland C4 function");
  createTag(*this, RBF_CPOLYNOMIAL_C6, once, SUBTAG_BASIS_FUNCTION, supportRadiusRBF, "Wendland C6 function");
  createTag(*this, RBF_CTPS_C2, once, SUBTAG_BASIS_FUNCTION, supportRadiusRBF, "Compact thin-plate-spline C2");
  createTag(*this, RBF_GAUSSIAN, once, SUBTAG_BASIS_FUNCTION, supportRadiusRBF, "Gaussian basis function accepting a support radius");

  auto attrSupportRadius = XMLAttribute<double>(ATTR_SUPPORT_RADIUS)
                               .setDocumentation("Support radius of each RBF basis function (global choice).");

  addAttributes(supportRadiusRBF, {attrSupportRadius});
  addSubtagsToParents(supportRadiusRBF, rbfIterativeTags);
  addSubtagsToParents(supportRadiusRBF, rbfDirectTags);

  // Now the tags using a shape parameter
  std::list<XMLTag> shapeParameterRBF;
  createTag(*this, RBF_MULTIQUADRICS, once, SUBTAG_BASIS_FUNCTION, shapeParameterRBF, "Multiquadrics");
  createTag(*this, RBF_INV_MULTIQUADRICS, once, SUBTAG_BASIS_FUNCTION, shapeParameterRBF, "Inverse multiquadrics");
  createTag(*this, RBF_GAUSSIAN, once, SUBTAG_BASIS_FUNCTION, shapeParameterRBF, "Gaussian basis function accepting a shape parameter");

  auto attrShapeParam = XMLAttribute<double>(ATTR_SHAPE_PARAM)
                            .setDocumentation("Specific shape parameter for RBF basis function.");

  addAttributes(shapeParameterRBF, {attrShapeParam});
  addSubtagsToParents(shapeParameterRBF, rbfIterativeTags);
  addSubtagsToParents(shapeParameterRBF, rbfDirectTags);

  // tags without an attribute
  std::list<XMLTag> attributelessRBFs;
  createTag(*this, RBF_TPS, once, SUBTAG_BASIS_FUNCTION, attributelessRBFs, "Thin-plate-splines");
  createTag(*this, RBF_VOLUME_SPLINES, once, SUBTAG_BASIS_FUNCTION, attributelessRBFs, "Volume splines");

  addSubtagsToParents(attributelessRBFs, rbfIterativeTags);
  addSubtagsToParents(attributelessRBFs, rbfDirectTags);

  // Add all tags to the mapping tag
  std::for_each(projectionTags.begin(), projectionTags.end(), [&parent](auto &s) { parent.addSubtag(s); });
  std::for_each(rbfIterativeTags.begin(), rbfIterativeTags.end(), [&parent](auto &s) { parent.addSubtag(s); });
  std::for_each(rbfDirectTags.begin(), rbfDirectTags.end(), [&parent](auto &s) { parent.addSubtag(s); });
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
    double        shapeParameter = std::numeric_limits<double>::quiet_NaN();
    double        supportRadius  = std::numeric_limits<double>::quiet_NaN();
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

    RBFParameter rbfParameter;
    // Check valid combinations for the Gaussian RBF input
    if (type == VALUE_RBF_GAUSSIAN || type == VALUE_RBF_INV_MULTIQUADRICS || type == VALUE_RBF_MULTIQUADRICS || type == VALUE_RBF_CTPS_C2 || type == VALUE_RBF_CPOLYNOMIAL_C0 || type == VALUE_RBF_CPOLYNOMIAL_C6 || type == VALUE_RBF_CPOLYNOMIAL_C2 || type == VALUE_RBF_CPOLYNOMIAL_C4) {
      const bool exactlyOneSet = (std::isfinite(supportRadius) && !std::isfinite(shapeParameter)) ||
                                 (std::isfinite(shapeParameter) && !std::isfinite(supportRadius));
      PRECICE_CHECK(exactlyOneSet, "The specified parameters for the Gaussian RBF mapping are invalid. Please specify either a \"shape-parameter\" or a \"support-radius\".");
      if (std::isfinite(shapeParameter)) {
        rbfParameter.type  = RBFParameter::Type::ShapeParameter;
        rbfParameter.value = shapeParameter;
      } else {
        rbfParameter.type  = RBFParameter::Type::SupportRadius;
        rbfParameter.value = supportRadius;
      }
    }
    ConfiguredMapping configuredMapping = createMapping(context,
                                                        dir, type, constraint,
                                                        fromMesh, toMesh,
                                                        rbfParameter, solverRtol,
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
    const RBFParameter &             rbfParameter,
    double                           solverRtol,
    bool                             xDead,
    bool                             yDead,
    bool                             zDead,
    bool                             useLU,
    Polynomial                       polynomial,
    Preallocation                    preallocation) const
{
  PRECICE_TRACE(direction, type, rbfParameter.value);
  using namespace mapping;
  ConfiguredMapping configuredMapping;
  mesh::PtrMesh     fromMesh(_meshConfig->getMesh(fromMeshName));
  mesh::PtrMesh     toMesh(_meshConfig->getMesh(toMeshName));
  PRECICE_CHECK(fromMesh.get() != nullptr,
                "Mesh \"{0}\" was not found while creating a mapping. "
                "Please correct the from=\"{0}\" attribute.",
                fromMeshName);
  PRECICE_CHECK(toMesh.get() != nullptr,
                "Mesh \"{0}\" was not found while creating a mapping. "
                "Please correct the to=\"{0}\" attribute.",
                toMeshName);
  configuredMapping.fromMesh = fromMesh;
  configuredMapping.toMesh   = toMesh;
  int dimensions             = fromMesh->getDimensions();

  if (direction == VALUE_WRITE) {
    configuredMapping.direction = WRITE;
  } else if (direction == VALUE_READ) {
    configuredMapping.direction = READ;
  } else {
    PRECICE_UNREACHABLE("Unknown mapping direction type \"{}\".", direction);
  }

  Mapping::Constraint constraintValue;
  if (constraint == VALUE_CONSERVATIVE) {
    constraintValue = Mapping::CONSERVATIVE;
  } else if (constraint == VALUE_CONSISTENT) {
    constraintValue = Mapping::CONSISTENT;
  } else if (constraint == VALUE_SCALED_CONSISTENT_SURFACE) {
    constraintValue = Mapping::SCALED_CONSISTENT_SURFACE;
  } else if (constraint == VALUE_SCALED_CONSISTENT_VOLUME) {
    constraintValue = Mapping::SCALED_CONSISTENT_VOLUME;
  } else {
    PRECICE_UNREACHABLE("Unknown mapping constraint \"{}\".", constraint);
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
  } else if (type == VALUE_LINEAR_CELL_INTERPOLATION) {
    configuredMapping.mapping = PtrMapping(
        new LinearCellInterpolationMapping(constraintValue, dimensions));
    configuredMapping.isRBF = false;
    return configuredMapping;
  } else if (type == VALUE_NEAREST_NEIGHBOR_GRADIENT) {

    // NNG is not applicable with the conservative constraint
    PRECICE_CHECK(constraintValue != Mapping::CONSERVATIVE,
                  "Nearest-neighbor-gradient mapping is not implemented using a \"conservative\" constraint. "
                  "Please select constraint=\" consistent\" or a different mapping method.");

    configuredMapping.mapping = PtrMapping(
        new NearestNeighborGradientMapping(constraintValue, dimensions));
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
          new RadialBasisFctMapping<ThinPlateSplines>(constraintValue, dimensions, ThinPlateSplines(), {{xDead, yDead, zDead}}, polynomial));
    } else if (type == VALUE_RBF_MULTIQUADRICS) {
      PRECICE_ASSERT(rbfParameter.type == RBFParameter::Type::ShapeParameter)
      configuredMapping.mapping = PtrMapping(
          new RadialBasisFctMapping<Multiquadrics>(
              constraintValue, dimensions, Multiquadrics(rbfParameter.value), {{xDead, yDead, zDead}}, polynomial));
    } else if (type == VALUE_RBF_INV_MULTIQUADRICS) {
      PRECICE_ASSERT(rbfParameter.type == RBFParameter::Type::ShapeParameter)
      configuredMapping.mapping = PtrMapping(
          new RadialBasisFctMapping<InverseMultiquadrics>(
              constraintValue, dimensions, InverseMultiquadrics(rbfParameter.value), {{xDead, yDead, zDead}}, polynomial));
    } else if (type == VALUE_RBF_VOLUME_SPLINES) {
      configuredMapping.mapping = PtrMapping(
          new RadialBasisFctMapping<VolumeSplines>(constraintValue, dimensions, VolumeSplines(), {{xDead, yDead, zDead}}, polynomial));
    } else if (type == VALUE_RBF_GAUSSIAN) {
      double shapeParameter = rbfParameter.value;
      if (rbfParameter.type == RBFParameter::Type::SupportRadius) {
        // Compute shape parameter from the support radius
        shapeParameter = std::sqrt(-std::log(Gaussian::cutoffThreshold)) / rbfParameter.value;
      }
      configuredMapping.mapping = PtrMapping(
          new RadialBasisFctMapping<Gaussian>(
              constraintValue, dimensions, Gaussian(shapeParameter), {{xDead, yDead, zDead}}, polynomial));
    } else if (type == VALUE_RBF_CTPS_C2) {
      PRECICE_ASSERT(rbfParameter.type == RBFParameter::Type::SupportRadius)
      configuredMapping.mapping = PtrMapping(
          new RadialBasisFctMapping<CompactThinPlateSplinesC2>(
              constraintValue, dimensions, CompactThinPlateSplinesC2(rbfParameter.value), {{xDead, yDead, zDead}}, polynomial));
    } else if (type == VALUE_RBF_CPOLYNOMIAL_C0) {
      PRECICE_ASSERT(rbfParameter.type == RBFParameter::Type::SupportRadius)
      configuredMapping.mapping = PtrMapping(
          new RadialBasisFctMapping<CompactPolynomialC0>(
              constraintValue, dimensions, CompactPolynomialC0(rbfParameter.value), {{xDead, yDead, zDead}}, polynomial));
    } else if (type == VALUE_RBF_CPOLYNOMIAL_C2) {
      PRECICE_ASSERT(rbfParameter.type == RBFParameter::Type::SupportRadius)
      configuredMapping.mapping = PtrMapping(
          new RadialBasisFctMapping<CompactPolynomialC2>(
              constraintValue, dimensions, CompactPolynomialC2(rbfParameter.value), {{xDead, yDead, zDead}}, polynomial));
    } else if (type == VALUE_RBF_CPOLYNOMIAL_C4) {
      PRECICE_ASSERT(rbfParameter.type == RBFParameter::Type::SupportRadius)
      configuredMapping.mapping = PtrMapping(
          new RadialBasisFctMapping<CompactPolynomialC4>(
              constraintValue, dimensions, CompactPolynomialC4(rbfParameter.value), {{xDead, yDead, zDead}}, polynomial));
    } else if (type == VALUE_RBF_CPOLYNOMIAL_C6) {
      PRECICE_ASSERT(rbfParameter.type == RBFParameter::Type::SupportRadius)
      configuredMapping.mapping = PtrMapping(
          new RadialBasisFctMapping<CompactPolynomialC6>(
              constraintValue, dimensions, CompactPolynomialC6(rbfParameter.value), {{xDead, yDead, zDead}}, polynomial));
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
                                                         {{xDead, yDead, zDead}}, solverRtol, polynomial, preallocation));
    } else if (type == VALUE_RBF_MULTIQUADRICS) {
      PRECICE_ASSERT(rbfParameter.type == RBFParameter::Type::ShapeParameter)
      configuredMapping.mapping = PtrMapping(
          new PetRadialBasisFctMapping<Multiquadrics>(constraintValue, dimensions, Multiquadrics(rbfParameter.value),
                                                      {{xDead, yDead, zDead}}, solverRtol, polynomial, preallocation));
    } else if (type == VALUE_RBF_INV_MULTIQUADRICS) {
      PRECICE_ASSERT(rbfParameter.type == RBFParameter::Type::ShapeParameter)
      configuredMapping.mapping = PtrMapping(
          new PetRadialBasisFctMapping<InverseMultiquadrics>(constraintValue, dimensions, InverseMultiquadrics(rbfParameter.value),
                                                             {{xDead, yDead, zDead}}, solverRtol, polynomial, preallocation));
    } else if (type == VALUE_RBF_VOLUME_SPLINES) {
      configuredMapping.mapping = PtrMapping(
          new PetRadialBasisFctMapping<VolumeSplines>(constraintValue, dimensions, VolumeSplines(),
                                                      {{xDead, yDead, zDead}}, solverRtol, polynomial, preallocation));
    } else if (type == VALUE_RBF_GAUSSIAN) {
      double shapeParameter = rbfParameter.value;
      if (rbfParameter.type == RBFParameter::Type::SupportRadius) {
        // Compute shape parameter from the support radius
        shapeParameter = std::sqrt(-std::log(Gaussian::cutoffThreshold)) / rbfParameter.value;
      }
      configuredMapping.mapping = PtrMapping(
          new PetRadialBasisFctMapping<Gaussian>(constraintValue, dimensions, Gaussian(shapeParameter),
                                                 {{xDead, yDead, zDead}}, solverRtol, polynomial, preallocation));
    } else if (type == VALUE_RBF_CTPS_C2) {
      PRECICE_ASSERT(rbfParameter.type == RBFParameter::Type::SupportRadius)
      configuredMapping.mapping = PtrMapping(
          new PetRadialBasisFctMapping<CompactThinPlateSplinesC2>(constraintValue, dimensions, CompactThinPlateSplinesC2(rbfParameter.value),
                                                                  {{xDead, yDead, zDead}}, solverRtol, polynomial, preallocation));
    } else if (type == VALUE_RBF_CPOLYNOMIAL_C0) {
      PRECICE_ASSERT(rbfParameter.type == RBFParameter::Type::SupportRadius)
      configuredMapping.mapping = PtrMapping(
          new PetRadialBasisFctMapping<CompactPolynomialC0>(constraintValue, dimensions, CompactPolynomialC0(rbfParameter.value),
                                                            {{xDead, yDead, zDead}}, solverRtol, polynomial, preallocation));
    } else if (type == VALUE_RBF_CPOLYNOMIAL_C2) {
      PRECICE_ASSERT(rbfParameter.type == RBFParameter::Type::SupportRadius)
      configuredMapping.mapping = PtrMapping(
          new PetRadialBasisFctMapping<CompactPolynomialC2>(constraintValue, dimensions, CompactPolynomialC2(rbfParameter.value),
                                                            {{xDead, yDead, zDead}}, solverRtol, polynomial, preallocation));
    } else if (type == VALUE_RBF_CPOLYNOMIAL_C4) {
      PRECICE_ASSERT(rbfParameter.type == RBFParameter::Type::SupportRadius)
      configuredMapping.mapping = PtrMapping(
          new PetRadialBasisFctMapping<CompactPolynomialC4>(constraintValue, dimensions, CompactPolynomialC4(rbfParameter.value),
                                                            {{xDead, yDead, zDead}}, solverRtol, polynomial, preallocation));
    } else if (type == VALUE_RBF_CPOLYNOMIAL_C6) {
      PRECICE_ASSERT(rbfParameter.type == RBFParameter::Type::SupportRadius)
      configuredMapping.mapping = PtrMapping(
          new PetRadialBasisFctMapping<CompactPolynomialC6>(constraintValue, dimensions, CompactPolynomialC6(rbfParameter.value),
                                                            {{xDead, yDead, zDead}}, solverRtol, polynomial, preallocation));
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
    PRECICE_CHECK(!sameMapping,
                  "There cannot be two mappings from mesh \"{}\" to mesh \"{}\". "
                  "Please remove one of the duplicated meshes. ",
                  mapping.fromMesh->getName(), mapping.toMesh->getName());
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

} // namespace precice::mapping
