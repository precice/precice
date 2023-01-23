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

// helper for the function below
template <typename>
inline constexpr bool always_false_v = false;

// Returns the value, if present
template <typename AttributeType>
AttributeType getAttributeIfPresent(xml::XMLTag &tag, const std::string &attributeName, AttributeType initializationValue)
{
  AttributeType value = initializationValue;
  if (tag.hasAttribute(attributeName)) {
    using T = std::decay_t<AttributeType>;
    if constexpr (std::is_same_v<T, double>) {
      value = tag.getDoubleAttributeValue(attributeName);
    } else if constexpr (std::is_same_v<T, bool>) {
      value = tag.getBooleanAttributeValue(attributeName);
    } else if constexpr (std::is_same_v<T, std::string>) {
      value = tag.getStringAttributeValue(attributeName);
    } else {
      static_assert(always_false_v<T>, "Attribute type is not supported.");
    }
  }
  return value;
}

enum struct RBFBackend {
  Eigen,
  PETSc
};

// Helper in order to resolve the template instantiations. Only the specializations are of interest
template <RBFBackend Backend, typename RBF>
struct BackendSelector {
  typedef RBF type;
};

template <typename RBF>
struct BackendSelector<RBFBackend::Eigen, RBF> {
  typedef mapping::RadialBasisFctMapping<RBF> type;
};

#ifndef PRECICE_NO_PETSC
template <typename RBF>
struct BackendSelector<RBFBackend::PETSc, RBF> {
  typedef mapping::PetRadialBasisFctMapping<RBF> type;
};
#endif

template <RBFBackend T, typename... Args>
PtrMapping instantiateRBFMapping(BasisFunctions functionType, mapping::Mapping::Constraint &constraint, int dimension, double supportRadius, double shapeParameter,
                                 Args &&... args)
{

  switch (functionType) {
  case BasisFunctions::WendlandC0: {
    return PtrMapping(new typename BackendSelector<T, CompactPolynomialC0>::type(constraint, dimension, mapping::CompactPolynomialC0(supportRadius), std::forward<Args>(args)...));
  }
  case BasisFunctions::WendlandC2: {
    return PtrMapping(new typename BackendSelector<T, CompactPolynomialC2>::type(constraint, dimension, mapping::CompactPolynomialC2(supportRadius), std::forward<Args>(args)...));
  }
  case BasisFunctions::WendlandC4: {
    return PtrMapping(new typename BackendSelector<T, CompactPolynomialC4>::type(constraint, dimension, mapping::CompactPolynomialC4(supportRadius), std::forward<Args>(args)...));
  }
  case BasisFunctions::WendlandC6: {
    return PtrMapping(new typename BackendSelector<T, CompactPolynomialC6>::type(constraint, dimension, mapping::CompactPolynomialC6(supportRadius), std::forward<Args>(args)...));
  }
  case BasisFunctions::CompactThinPlateSplinesC2: {
    return PtrMapping(new typename BackendSelector<T, CompactThinPlateSplinesC2>::type(constraint, dimension, mapping::CompactThinPlateSplinesC2(supportRadius), std::forward<Args>(args)...));
  }
  case BasisFunctions::ThinPlateSplines: {
    return PtrMapping(new typename BackendSelector<T, ThinPlateSplines>::type(constraint, dimension, mapping::ThinPlateSplines(), std::forward<Args>(args)...));
  }
  case BasisFunctions::Multiquadrics: {
    return PtrMapping(new typename BackendSelector<T, Multiquadrics>::type(constraint, dimension, mapping::Multiquadrics(shapeParameter), std::forward<Args>(args)...));
  }
  case BasisFunctions::InverseMultiquadrics: {
    return PtrMapping(new typename BackendSelector<T, InverseMultiquadrics>::type(constraint, dimension, mapping::InverseMultiquadrics(shapeParameter), std::forward<Args>(args)...));
  }
  case BasisFunctions::Gaussian: {
    return PtrMapping(new typename BackendSelector<T, Gaussian>::type(constraint, dimension, mapping::Gaussian(shapeParameter), std::forward<Args>(args)...));
  }
  default:
    PRECICE_UNREACHABLE("Unkown basis function.");
    return nullptr;
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
  createTag(*this, RBF_GAUSSIAN_SUPPORT, once, SUBTAG_BASIS_FUNCTION, supportRadiusRBF, "Gaussian basis function accepting a support radius");

  auto attrSupportRadius = XMLAttribute<double>(ATTR_SUPPORT_RADIUS)
                               .setDocumentation("Support radius of each RBF basis function (global choice).");

  addAttributes(supportRadiusRBF, {attrSupportRadius});
  addSubtagsToParents(supportRadiusRBF, rbfIterativeTags);
  addSubtagsToParents(supportRadiusRBF, rbfDirectTags);

  // Now the tags using a shape parameter
  std::list<XMLTag> shapeParameterRBF;
  createTag(*this, RBF_MULTIQUADRICS, once, SUBTAG_BASIS_FUNCTION, shapeParameterRBF, "Multiquadrics");
  createTag(*this, RBF_INV_MULTIQUADRICS, once, SUBTAG_BASIS_FUNCTION, shapeParameterRBF, "Inverse multiquadrics");
  createTag(*this, RBF_GAUSSIAN_SHAPE, once, SUBTAG_BASIS_FUNCTION, shapeParameterRBF, "Gaussian basis function accepting a shape parameter");

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
    // Mandatory tags
    std::string dir        = tag.getStringAttributeValue(ATTR_DIRECTION);
    std::string fromMesh   = tag.getStringAttributeValue(ATTR_FROM);
    std::string toMesh     = tag.getStringAttributeValue(ATTR_TO);
    std::string type       = tag.getName();
    std::string constraint = tag.getStringAttributeValue(ATTR_CONSTRAINT);

    // optional tags
    bool        xDead         = getAttributeIfPresent(tag, ATTR_X_DEAD, false);
    bool        yDead         = getAttributeIfPresent(tag, ATTR_Y_DEAD, false);
    bool        zDead         = getAttributeIfPresent(tag, ATTR_Z_DEAD, false);
    double      solverRtol    = getAttributeIfPresent(tag, ATTR_SOLVER_RTOL, 1e-9);
    std::string strPolynomial = getAttributeIfPresent(tag, ATTR_POLYNOMIAL, POLYNOMIAL_SEPARATE);
    std::string strPrealloc   = getAttributeIfPresent(tag, ATTR_PREALLOCATION, PREALLOCATION_TREE);

    // Convert raw string into enum types as the constructors take enums
    if (constraint == CONSTRAINT_CONSERVATIVE) {
      constraintValue = Mapping::CONSERVATIVE;
    } else if (constraint == CONSTRAINT_CONSISTENT) {
      constraintValue = Mapping::CONSISTENT;
    } else if (constraint == CONSTRAINT_SCALED_CONSISTENT_SURFACE) {
      constraintValue = Mapping::SCALED_CONSISTENT_SURFACE;
    } else if (constraint == CONSTRAINT_SCALED_CONSISTENT_VOLUME) {
      constraintValue = Mapping::SCALED_CONSISTENT_VOLUME;
    } else {
      PRECICE_UNREACHABLE("Unknown mapping constraint \"{}\".", constraint);
    }

    ConfiguredMapping configuredMapping = createMapping(dir, type, fromMesh, toMesh);

    _rbfConfig = configureRBFMapping(type, strPolynomial, strPrealloc, xDead, yDead, zDead, solverRtol);

    checkDuplicates(configuredMapping);
    _mappings.push_back(configuredMapping);
  } else if (tag.getNamespace() == SUBTAG_BASIS_FUNCTION) {

    PRECICE_ASSERT(!_mappings.empty());
    // We can only set one subtag
    PRECICE_CHECK(_mappings.back().mapping == nullptr, "More than one basis-function was defined for the.");

    std::string basisFctName   = tag.getName();
    double      supportRadius  = getAttributeIfPresent(tag, ATTR_SUPPORT_RADIUS, 0.);
    double      shapeParameter = getAttributeIfPresent(tag, ATTR_SHAPE_PARAM, 0.);

    ConfiguredMapping &mapping = _mappings.back();

    BasisFunctions basisFunction{};
    if (basisFctName == RBF_TPS)
      basisFunction = BasisFunctions::ThinPlateSplines;
    else if (basisFctName == RBF_MULTIQUADRICS)
      basisFunction = BasisFunctions::Multiquadrics;
    else if (basisFctName == RBF_INV_MULTIQUADRICS)
      basisFunction = BasisFunctions::InverseMultiquadrics;
    else if (basisFctName == RBF_VOLUME_SPLINES)
      basisFunction = BasisFunctions::VolumeSplines;
    else if (basisFctName == RBF_GAUSSIAN_SHAPE || basisFctName == RBF_GAUSSIAN_SUPPORT)
      basisFunction = BasisFunctions::Gaussian;
    else if (basisFctName == RBF_CTPS_C2)
      basisFunction = BasisFunctions::CompactThinPlateSplinesC2;
    else if (basisFctName == RBF_CPOLYNOMIAL_C0)
      basisFunction = BasisFunctions::WendlandC0;
    else if (basisFctName == RBF_CPOLYNOMIAL_C2)
      basisFunction = BasisFunctions::WendlandC2;
    else if (basisFctName == RBF_CPOLYNOMIAL_C4)
      basisFunction = BasisFunctions::WendlandC4;
    else if (basisFctName == RBF_CPOLYNOMIAL_C6)
      basisFunction = BasisFunctions::WendlandC6;
    else
      PRECICE_UNREACHABLE("Unknown basis function \"{}\".", basisFctName);

    // The Gaussian RBF is always treated as a shape-parameter RBF. Hence, we have to convert the support radius, if necessary
    if (basisFunction == BasisFunctions::Gaussian && supportRadius > 0 && shapeParameter == 0) {
      shapeParameter = std::sqrt(-std::log(Gaussian::cutoffThreshold)) / supportRadius;
    }

    if (_rbfConfig.solver == RBFConfiguration::SystemSolver::GlobalDirect) {
      mapping.mapping = instantiateRBFMapping<RBFBackend::Eigen>(basisFunction, constraintValue, mapping.fromMesh->getDimensions(), supportRadius, shapeParameter, _rbfConfig.deadAxis, _rbfConfig.polynomial);
    } else if (_rbfConfig.solver == RBFConfiguration::SystemSolver::GlobalIterative) {
#ifndef PRECICE_NO_PETSC
      // for petsc initialization
      int   argc = 1;
      char *arg  = new char[8];
      strcpy(arg, "precice");
      char **argv = &arg;
      utils::Petsc::initialize(&argc, &argv, utils::Parallel::current()->comm);
      delete[] arg;

      mapping.mapping = instantiateRBFMapping<RBFBackend::PETSc>(basisFunction, constraintValue, mapping.fromMesh->getDimensions(), supportRadius, shapeParameter, _rbfConfig.deadAxis, _rbfConfig.solverRtol, _rbfConfig.polynomial, _rbfConfig.preallocation);
#else
      PRECICE_CHECK(false, "The global-iterative RBF solver requires a preCICE build with PETSc enabled.");
#endif
    } else {
      PRECICE_UNREACHABLE("Unknown RBF solver.");
    }
  }
}

MappingConfiguration::RBFConfiguration MappingConfiguration::configureRBFMapping(const std::string &type,
                                                                                 const std::string &polynomial,
                                                                                 const std::string &preallocation,
                                                                                 bool xDead, bool yDead, bool zDead,
                                                                                 double solverRtol) const
{
  RBFConfiguration rbfConfig;

  if (type == TYPE_RBF_GLOBAL_ITERATIVE)
    rbfConfig.solver = RBFConfiguration::SystemSolver::GlobalIterative;
  else if (type == TYPE_RBF_GLOBAL_DIRECT)
    rbfConfig.solver = RBFConfiguration::SystemSolver::GlobalDirect;

  if (polynomial == POLYNOMIAL_SEPARATE)
    rbfConfig.polynomial = Polynomial::SEPARATE;
  else if (polynomial == POLYNOMIAL_ON)
    rbfConfig.polynomial = Polynomial::ON;
  else if (polynomial == POLYNOMIAL_OFF)
    rbfConfig.polynomial = Polynomial::OFF;

  if (preallocation == PREALLOCATION_ESTIMATE)
    rbfConfig.preallocation = Preallocation::ESTIMATE;
  else if (preallocation == PREALLOCATION_COMPUTE)
    rbfConfig.preallocation = Preallocation::COMPUTE;
  else if (preallocation == PREALLOCATION_SAVE)
    rbfConfig.preallocation = Preallocation::SAVE;
  else if (preallocation == PREALLOCATION_TREE)
    rbfConfig.preallocation = Preallocation::TREE;
  else if (preallocation == PREALLOCATION_OFF)
    rbfConfig.preallocation = Preallocation::OFF;

  rbfConfig.deadAxis   = {{xDead, yDead, zDead}};
  rbfConfig.solverRtol = solverRtol;

  return rbfConfig;
}

MappingConfiguration::ConfiguredMapping MappingConfiguration::createMapping(
    const std::string &direction,
    const std::string &type,
    const std::string &fromMeshName,
    const std::string &toMeshName) const
{
  PRECICE_TRACE(direction, type);

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

  if (direction == DIRECTION_WRITE) {
    configuredMapping.direction = WRITE;
  } else if (direction == DIRECTION_READ) {
    configuredMapping.direction = READ;
  } else {
    PRECICE_UNREACHABLE("Unknown mapping direction type \"{}\".", direction);
  }

  // Create all projection based mappings
  if (type == TYPE_NEAREST_NEIGHBOR) {
    configuredMapping.mapping = PtrMapping(new NearestNeighborMapping(constraintValue, fromMesh->getDimensions()));
  } else if (type == TYPE_NEAREST_PROJECTION) {
    configuredMapping.mapping = PtrMapping(new NearestProjectionMapping(constraintValue, fromMesh->getDimensions()));
  } else if (type == TYPE_LINEAR_CELL_INTERPOLATION) {
    configuredMapping.mapping = PtrMapping(new LinearCellInterpolationMapping(constraintValue, fromMesh->getDimensions()));
  } else if (type == TYPE_NEAREST_NEIGHBOR_GRADIENT) {

    // NNG is not applicable with the conservative constraint
    PRECICE_CHECK(constraintValue != Mapping::CONSERVATIVE,
                  "Nearest-neighbor-gradient mapping is not implemented using a \"conservative\" constraint. "
                  "Please select constraint=\" consistent\" or a different mapping method.");

    configuredMapping.mapping = PtrMapping(new NearestNeighborGradientMapping(constraintValue, fromMesh->getDimensions()));
  } else {
    // We need knowledge about the basis function in order to instantiate the rbf related mapping
    PRECICE_ASSERT(requiresBasisFunction(type));
    configuredMapping.mapping = nullptr;
  }

  configuredMapping.requiresBasisFunction = requiresBasisFunction(type);

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
  if (requiresBasisFunction(tag.getName())) {
    PRECICE_CHECK(_mappings.back().mapping != nullptr, "No basis-function was defined for the \"{}\" mapping from mesh \"{}\" to mesh \"{}\".", tag.getName(), _mappings.back().fromMesh->getName(), _mappings.back().toMesh->getName());
  }
  PRECICE_ASSERT(_mappings.back().mapping != nullptr);
}

const std::vector<MappingConfiguration::ConfiguredMapping> &MappingConfiguration::mappings()
{
  return _mappings;
}

bool MappingConfiguration::requiresBasisFunction(const std::string &mappingType) const
{
  return mappingType == TYPE_RBF_GLOBAL_DIRECT || mappingType == TYPE_RBF_GLOBAL_ITERATIVE;
}
} // namespace precice::mapping
