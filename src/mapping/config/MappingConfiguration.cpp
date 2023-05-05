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
#include "mapping/AxialGeoMultiscaleMapping.hpp"
#include "mapping/LinearCellInterpolationMapping.hpp"
#include "mapping/Mapping.hpp"
#include "mapping/NearestNeighborGradientMapping.hpp"
#include "mapping/NearestNeighborMapping.hpp"
#include "mapping/NearestProjectionMapping.hpp"
#include "mapping/PartitionOfUnityMapping.hpp"
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

// given a list of subtags and parent tags, this function adds all subtags to all
// parent tags
void addSubtagsToParents(std::list<xml::XMLTag> &subtags,
                         std::list<xml::XMLTag> &parents)
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

  auto attrRadius = XMLAttribute<double>("radius")
                        .setDocumentation("Radius for 1D participants in a geometric multiscale mapping.");
  auto attrMultiscaleType = XMLAttribute<std::string>("type")
                                .setDocumentation("Type of a geometric multiscale mapping (spread or collect).")
                                .setOptions({"spread", "collect"});

  XMLTag::Occurrence occ = XMLTag::OCCUR_ARBITRARY;
  std::list<XMLTag>  tags;
  {
    XMLTag tag(*this, VALUE_RBF_TPS, occ, TAG);
    tag.setDocumentation("Global radial-basis-function mapping based on the thin plate splines.");
    tags.push_back(tag);
  }
}

// this function uses std::variant in order to add attributes of any type (double, string, bool)
// to all tags in the list of tags \p storage.
using variant_t = std::variant<xml::XMLAttribute<double>, xml::XMLAttribute<std::string>, xml::XMLAttribute<bool>, xml::XMLAttribute<int>>;
template <typename TagStorage>
void addAttributes(TagStorage &storage, const std::vector<variant_t> &attributes)
{
  for (auto &s : storage) {
    for (auto &a : attributes)
      std::visit([&s](auto &&arg) { s.addAttribute(arg); }, a);
  }
}

// Enum required for the RBF instantiations
enum struct RBFBackend {
  Eigen,
  PETSc,
  PUM
};

// Helper in order to resolve the template instantiations.
// Only the template specializations are of interest
template <RBFBackend Backend, typename RBF>
struct BackendSelector {
  typedef RBF type;
};

// Specialization for the RBF Eigen backend
template <typename RBF>
struct BackendSelector<RBFBackend::Eigen, RBF> {
  typedef mapping::RadialBasisFctMapping<RBF> type;
};

// Specialization for the PETSc RBF backend
#ifndef PRECICE_NO_PETSC
template <typename RBF>
struct BackendSelector<RBFBackend::PETSc, RBF> {
  typedef mapping::PetRadialBasisFctMapping<RBF> type;
};
#endif

// Specialization for the RBF PUM backend
template <typename RBF>
struct BackendSelector<RBFBackend::PUM, RBF> {
  typedef mapping::PartitionOfUnityMapping<RBF> type;
};

// Variant holding all available RBF classes
using rbf_variant_t = std::variant<CompactPolynomialC0, CompactPolynomialC2, CompactPolynomialC4, CompactPolynomialC6, CompactThinPlateSplinesC2, ThinPlateSplines, VolumeSplines, Multiquadrics, InverseMultiquadrics, Gaussian>;

// The actual instantiation of the mapping class, which is called by the visitor \ref getRBFMapping
template <RBFBackend T, typename RADIAL_BASIS_FUNCTION_T, typename... Args>
PtrMapping instantiateRBFMapping(mapping::Mapping::Constraint &constraint, int dimension, RADIAL_BASIS_FUNCTION_T function,
                                 Args &&... args)
{
  return PtrMapping(new typename BackendSelector<T, RADIAL_BASIS_FUNCTION_T>::type(constraint, dimension, function, std::forward<Args>(args)...));
}

// Constructs the RBF function based on the functionType
rbf_variant_t constructRBF(BasisFunction functionType, double supportRadius, double shapeParameter)
{
  switch (functionType) {
  case BasisFunction::WendlandC0: {
    return mapping::CompactPolynomialC0(supportRadius);
  }
  case BasisFunction::WendlandC2: {
    return mapping::CompactPolynomialC2(supportRadius);
  }
  case BasisFunction::WendlandC4: {
    return mapping::CompactPolynomialC4(supportRadius);
  }
  case BasisFunction::WendlandC6: {
    return mapping::CompactPolynomialC6(supportRadius);
  }
  case BasisFunction::CompactThinPlateSplinesC2: {
    return mapping::CompactThinPlateSplinesC2(supportRadius);
  }
  case BasisFunction::ThinPlateSplines: {
    return mapping::ThinPlateSplines();
  }
  case BasisFunction::VolumeSplines: {
    return mapping::VolumeSplines();
  }
  case BasisFunction::Multiquadrics: {
    return mapping::Multiquadrics(shapeParameter);
  }
  case BasisFunction::InverseMultiquadrics: {
    return mapping::InverseMultiquadrics(shapeParameter);
  }
  case BasisFunction::Gaussian: {
    return mapping::Gaussian(shapeParameter);
  }
  default:
    PRECICE_UNREACHABLE("No instantiation was found for the selected basis function.");
  }
  {
    XMLTag tag(*this, VALUE_AXIAL_GEOMETRIC_MULTISCALE, occ, TAG);
    tag.setDocumentation("Axial geometric multiscale mapping.");
    tag.addAttribute(attrRadius);
    tag.addAttribute(attrMultiscaleType);
    tags.push_back(tag);
  }

// The actual instantion helper, which avoids enumerating all mapping implementations (more will come) with all RBF kernels
// The first three arguments of the constructor are prescribed: constraint, dimension and the RBF function object, all other
// constructor arguments are just forwareded. The first argument (BasisFunction) indicates then the actual instantiation to return.
template <RBFBackend T, typename... Args>
PtrMapping getRBFMapping(BasisFunction functionType, mapping::Mapping::Constraint &constraint, int dimension, double supportRadius, double shapeParameter,
                         Args &&... args)
{
  // First, construct the RBF function
  auto functionVariant = constructRBF(functionType, supportRadius, shapeParameter);
  // ... and instantiate the corresponding RBF mapping class
  return std::visit([&](auto &&func) { return instantiateRBFMapping<T>(constraint, dimension, func, std::forward<Args>(args)...); }, functionVariant);
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
  std::list<XMLTag>  projectionTags{
      XMLTag{*this, TYPE_NEAREST_NEIGHBOR, occ, TAG}.setDocumentation("Nearest-neighbour mapping which uses a rstar-spacial index tree to index meshes and run nearest-neighbour queries."),
      XMLTag{*this, TYPE_NEAREST_PROJECTION, occ, TAG}.setDocumentation("Nearest-projection mapping which uses a rstar-spacial index tree to index meshes and locate the nearest projections."),
      XMLTag{*this, TYPE_NEAREST_NEIGHBOR_GRADIENT, occ, TAG}.setDocumentation("Nearest-neighbor-gradient mapping which uses nearest-neighbor mapping with an additional linear approximation using gradient data."),
      XMLTag{*this, TYPE_LINEAR_CELL_INTERPOLATION, occ, TAG}.setDocumentation("Linear cell interpolation mapping which uses a rstar-spacial index tree to index meshes and locate the nearest cell. Only supports 2D meshes.")};
  std::list<XMLTag> rbfDirectTags{
      XMLTag{*this, TYPE_RBF_GLOBAL_DIRECT, occ, TAG}.setDocumentation("Radial-basis-function mapping using a direct solver with a gather-scatter parallelism.")};
  std::list<XMLTag> rbfIterativeTags{
      XMLTag{*this, TYPE_RBF_GLOBAL_ITERATIVE, occ, TAG}.setDocumentation("Radial-basis-function mapping using an iterative solver with a distributed parallelism.")};
  std::list<XMLTag> pumDirectTags{
      XMLTag{*this, TYPE_RBF_PUM_DIRECT, occ, TAG}.setDocumentation("Radial-basis-function mapping using a partition of unity method, which supports a distributed parallelism.")};
  std::list<XMLTag> rbfAliasTag{
      XMLTag{*this, TYPE_RBF_ALIAS, occ, TAG}.setDocumentation("Alias tag, which auto-selects a radial-basis-function mapping depending on the simulation parameter,")};

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
  auto attrPumPolynomial = makeXMLAttribute(ATTR_POLYNOMIAL, POLYNOMIAL_SEPARATE)
                               .setDocumentation("Toggles use a local (per cluster) polynomial")
                               .setOptions({POLYNOMIAL_OFF, POLYNOMIAL_SEPARATE});

  auto attrSolverRtol = makeXMLAttribute(ATTR_SOLVER_RTOL, 1e-9)
                            .setDocumentation("Solver relative tolerance for convergence");
  auto attrPreallocation = makeXMLAttribute(ATTR_PREALLOCATION, PREALLOCATION_TREE)
                               .setDocumentation("Sets kind of preallocation for PETSc RBF implementation")
                               .setOptions({PREALLOCATION_ESTIMATE, PREALLOCATION_COMPUTE, PREALLOCATION_OFF, PREALLOCATION_SAVE, PREALLOCATION_TREE});

  auto verticesPerCluster = XMLAttribute<int>(ATTR_VERTICES_PER_CLUSTER, 100)
                                .setDocumentation("Average number of vertices per cluster (partition) applied in the rbf partition of unity method.");
  auto relativeOverlap = makeXMLAttribute(ATTR_RELATIVE_OVERLAP, 0.3)
                             .setDocumentation("Value between 0 and 1 indicating the relative overlap between clusters. A value of 0.3 is usually a good trade-off between accuracy and efficiency.");
  auto projectToInput = XMLAttribute<bool>(ATTR_PROJECT_TO_INPUT, true)
                            .setDocumentation("If enabled, places the cluster centers at the closest vertex of the input mesh. Should be enabled in case of non-uniform point distributions such as for shell structures.");

  // Add the relevant attributes to the relevant tags
  addAttributes(projectionTags, {attrFromMesh, attrToMesh, attrDirection, attrConstraint});
  addAttributes(rbfDirectTags, {attrFromMesh, attrToMesh, attrDirection, attrConstraint, attrPolynomial, attrXDead, attrYDead, attrZDead});
  addAttributes(rbfIterativeTags, {attrFromMesh, attrToMesh, attrDirection, attrConstraint, attrPolynomial, attrXDead, attrYDead, attrZDead, attrSolverRtol, attrPreallocation});
  addAttributes(pumDirectTags, {attrFromMesh, attrToMesh, attrDirection, attrConstraint, attrPumPolynomial, attrXDead, attrYDead, attrZDead, verticesPerCluster, relativeOverlap, projectToInput});
  addAttributes(rbfAliasTag, {attrFromMesh, attrToMesh, attrDirection, attrConstraint, attrXDead, attrYDead, attrZDead});

  // Now we take care of the subtag basis function
  // First, we have the tags using a support radius
  XMLTag::Occurrence once = XMLTag::OCCUR_NOT_OR_ONCE;
  std::list<XMLTag>  supportRadiusRBF{
      XMLTag{*this, RBF_CPOLYNOMIAL_C0, once, SUBTAG_BASIS_FUNCTION}.setDocumentation("Wendland C0 function"),
      XMLTag{*this, RBF_CPOLYNOMIAL_C2, once, SUBTAG_BASIS_FUNCTION}.setDocumentation("Wendland C2 function"),
      XMLTag{*this, RBF_CPOLYNOMIAL_C4, once, SUBTAG_BASIS_FUNCTION}.setDocumentation("Wendland C4 function"),
      XMLTag{*this, RBF_CPOLYNOMIAL_C6, once, SUBTAG_BASIS_FUNCTION}.setDocumentation("Wendland C6 function"),
      XMLTag{*this, RBF_CTPS_C2, once, SUBTAG_BASIS_FUNCTION}.setDocumentation("Compact thin-plate-spline C2")};

  auto attrSupportRadius = XMLAttribute<double>(ATTR_SUPPORT_RADIUS)
                               .setDocumentation("Support radius of each RBF basis function (global choice).");

  addAttributes(supportRadiusRBF, {attrSupportRadius});
  addSubtagsToParents(supportRadiusRBF, rbfIterativeTags);
  addSubtagsToParents(supportRadiusRBF, rbfDirectTags);
  addSubtagsToParents(supportRadiusRBF, pumDirectTags);
  addSubtagsToParents(supportRadiusRBF, rbfAliasTag);

  // Now the tags using a shape parameter
  std::list<XMLTag> shapeParameterRBF{
      XMLTag{*this, RBF_MULTIQUADRICS, once, SUBTAG_BASIS_FUNCTION}.setDocumentation("Multiquadrics"),
      XMLTag{*this, RBF_INV_MULTIQUADRICS, once, SUBTAG_BASIS_FUNCTION}.setDocumentation("Inverse multiquadrics")};

  auto attrShapeParam = XMLAttribute<double>(ATTR_SHAPE_PARAM)
                            .setDocumentation("Specific shape parameter for RBF basis function.");

  addAttributes(shapeParameterRBF, {attrShapeParam});
  addSubtagsToParents(shapeParameterRBF, rbfIterativeTags);
  addSubtagsToParents(shapeParameterRBF, rbfDirectTags);
  addSubtagsToParents(shapeParameterRBF, pumDirectTags);
  addSubtagsToParents(shapeParameterRBF, rbfAliasTag);

  // For the Gaussian, we need default values as the user can pass a support radius or a shape parameter
  std::list<XMLTag> GaussRBF{
      XMLTag{*this, RBF_GAUSSIAN, once, SUBTAG_BASIS_FUNCTION}.setDocumentation("Gaussian basis function accepting a support radius or a shape parameter.")};
  attrShapeParam.setDefaultValue(std::numeric_limits<double>::quiet_NaN());
  attrSupportRadius.setDefaultValue(std::numeric_limits<double>::quiet_NaN());
  addAttributes(GaussRBF, {attrShapeParam, attrSupportRadius});
  addSubtagsToParents(GaussRBF, rbfIterativeTags);
  addSubtagsToParents(GaussRBF, rbfDirectTags);
  addSubtagsToParents(GaussRBF, pumDirectTags);
  addSubtagsToParents(GaussRBF, rbfAliasTag);

  // tags without an attribute
  std::list<XMLTag> attributelessRBFs{
      XMLTag{*this, RBF_TPS, once, SUBTAG_BASIS_FUNCTION}.setDocumentation("Thin-plate-splines"),
      XMLTag{*this, RBF_VOLUME_SPLINES, once, SUBTAG_BASIS_FUNCTION}.setDocumentation("Volume splines")};

  addSubtagsToParents(attributelessRBFs, rbfIterativeTags);
  addSubtagsToParents(attributelessRBFs, rbfDirectTags);
  addSubtagsToParents(attributelessRBFs, pumDirectTags);
  addSubtagsToParents(attributelessRBFs, rbfAliasTag);

  // Add all tags to the mapping tag
  parent.addSubtags(projectionTags);
  parent.addSubtags(rbfIterativeTags);
  parent.addSubtags(rbfDirectTags);
  parent.addSubtags(pumDirectTags);
  parent.addSubtags(rbfAliasTag);
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
    double        shapeParameter = std::numeric_limits<double>::quiet_NaN();
    double        supportRadius  = std::numeric_limits<double>::quiet_NaN();
    double        solverRtol     = 1e-9;
    bool          xDead = false, yDead = false, zDead = false;
    bool          useLU          = false;
    Polynomial    polynomial     = Polynomial::ON;
    Preallocation preallocation  = Preallocation::TREE;
    double        radius         = 0.0;
    std::string   multiscaleType = "undefined";

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
    if (tag.hasAttribute("radius")) {
      radius = tag.getDoubleAttributeValue("radius");
    }
    if (tag.hasAttribute("type")) {
      multiscaleType = tag.getStringAttributeValue("type");
    }

    ConfiguredMapping configuredMapping = createMapping(dir, type, fromMesh, toMesh);

    _rbfConfig = configureRBFMapping(type, strPolynomial, strPrealloc, xDead, yDead, zDead, solverRtol, verticesPerCluster, relativeOverlap, projectToInput);

    checkDuplicates(configuredMapping);
    _mappings.push_back(configuredMapping);
  } else if (tag.getNamespace() == SUBTAG_BASIS_FUNCTION) {

    PRECICE_ASSERT(!_mappings.empty());
    // We can only set one subtag
    PRECICE_CHECK(_mappings.back().mapping == nullptr, "More than one basis-function was defined for the mapping "
                                                       "from mesh \"{}\" to mesh \"{}\".",
                  _mappings.back().fromMesh->getName(), _mappings.back().toMesh->getName());

    std::string basisFctName   = tag.getName();
    double      supportRadius  = tag.getDoubleAttributeValue(ATTR_SUPPORT_RADIUS, 0.);
    double      shapeParameter = tag.getDoubleAttributeValue(ATTR_SHAPE_PARAM, 0.);

    ConfiguredMapping &mapping = _mappings.back();

    BasisFunction basisFunction = parseBasisFunctions(basisFctName);

    // The Gaussian RBF is always treated as a shape-parameter RBF. Hence, we have to convert the support radius, if necessary
    if (basisFunction == BasisFunction::Gaussian) {
      const bool exactlyOneSet = (std::isfinite(supportRadius) && !std::isfinite(shapeParameter)) ||
                                 (std::isfinite(shapeParameter) && !std::isfinite(supportRadius));
      PRECICE_CHECK(exactlyOneSet, "The specified parameters for the Gaussian RBF mapping are invalid. Please specify either a \"shape-parameter\" or a \"support-radius\".");

      if (std::isfinite(supportRadius) && !std::isfinite(shapeParameter)) {
        shapeParameter = std::sqrt(-std::log(Gaussian::cutoffThreshold)) / supportRadius;
      }
    }
    ConfiguredMapping configuredMapping = createMapping(context,
                                                        dir, type, constraint,
                                                        fromMesh, toMesh, timing,
                                                        rbfParameter, solverRtol,
                                                        xDead, yDead, zDead,
                                                        useLU,
                                                        polynomial, preallocation,
                                                        radius, multiscaleType);
    checkDuplicates(configuredMapping);
    _mappings.push_back(configuredMapping);
  }
}

MappingConfiguration::RBFConfiguration MappingConfiguration::configureRBFMapping(const std::string &type,
                                                                                 const std::string &polynomial,
                                                                                 const std::string &preallocation,
                                                                                 bool xDead, bool yDead, bool zDead,
                                                                                 double solverRtol,
                                                                                 double verticesPerCluster,
                                                                                 double relativeOverlap,
                                                                                 bool   projectToInput) const
{
  RBFConfiguration rbfConfig;

  if (type == TYPE_RBF_GLOBAL_ITERATIVE)
    rbfConfig.solver = RBFConfiguration::SystemSolver::GlobalIterative;
  else if (type == TYPE_RBF_GLOBAL_DIRECT)
    rbfConfig.solver = RBFConfiguration::SystemSolver::GlobalDirect;
  else if (type == TYPE_RBF_PUM_DIRECT)
    rbfConfig.solver = RBFConfiguration::SystemSolver::PUMDirect;
  else {
    // Rather simple auto-selection (for now), consisting of the PUM Eigen backend
    rbfConfig.solver = RBFConfiguration::SystemSolver::PUMDirect;

    // A more sophisticated criterion here could take the globalNumberOfVertices into account
    // (the mesh pointer is stored in the configuredMapping anyway), but this quantity is not yet computed
    // during the configuration time.
  }

  if (polynomial == POLYNOMIAL_SEPARATE)
    rbfConfig.polynomial = Polynomial::SEPARATE;
  else if (polynomial == POLYNOMIAL_ON)
    rbfConfig.polynomial = Polynomial::ON;
  else if (polynomial == POLYNOMIAL_OFF)
    rbfConfig.polynomial = Polynomial::OFF;
  else
    PRECICE_UNREACHABLE("Unknown polynomial configuration.");

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
  else
    PRECICE_UNREACHABLE("Unknown preallocation configuration");

  rbfConfig.deadAxis   = {{xDead, yDead, zDead}};
  rbfConfig.solverRtol = solverRtol;

  rbfConfig.verticesPerCluster = verticesPerCluster;
  rbfConfig.relativeOverlap    = relativeOverlap;
  rbfConfig.projectToInput     = projectToInput;

  return rbfConfig;
}

MappingConfiguration::ConfiguredMapping MappingConfiguration::createMapping(
    const xml::ConfigurationContext &context,
    const std::string &              direction,
    const std::string &              type,
    const std::string &              constraint,
    const std::string &              fromMeshName,
    const std::string &              toMeshName,
    Timing                           timing,
    const RBFParameter &             rbfParameter,
    double                           solverRtol,
    bool                             xDead,
    bool                             yDead,
    bool                             zDead,
    bool                             useLU,
    Polynomial                       polynomial,
    Preallocation                    preallocation,
    double                           radius,
    const std::string &              multiscaleType) const
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
    } else if (type == VALUE_RBF_CPOLYNOMIAL_C6) {
      PRECICE_ASSERT(rbfParameter.type == RBFParameter::Type::SupportRadius)
      configuredMapping.mapping = PtrMapping(
          new RadialBasisFctMapping<CompactPolynomialC6>(
              constraintValue, dimensions, CompactPolynomialC6(rbfParameter.value), {{xDead, yDead, zDead}}, polynomial));
    } else if (type == VALUE_AXIAL_GEOMETRIC_MULTISCALE) {

      AxialGeoMultiscaleMapping::MultiscaleType multiscaleTypeValue;
      if (multiscaleType == "spread") {
        multiscaleTypeValue = AxialGeoMultiscaleMapping::SPREAD;
      } else if (multiscaleType == "collect") {
        multiscaleTypeValue = AxialGeoMultiscaleMapping::COLLECT;
      } else {
        PRECICE_ERROR("Unknown geometric multiscale type \"{}\". Known types are \"spread\" and \"collect\".", multiscaleTypeValue);
      }
      configuredMapping.mapping = PtrMapping(
          new AxialGeoMultiscaleMapping(constraintValue, dimensions, multiscaleTypeValue, radius));
      configuredMapping.isRBF = false;
    } else {
      std::cout << type;
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
    } else if (type == VALUE_RBF_CPOLYNOMIAL_C6) {
      PRECICE_ASSERT(rbfParameter.type == RBFParameter::Type::SupportRadius)
      configuredMapping.mapping = PtrMapping(
          new PetRadialBasisFctMapping<CompactPolynomialC6>(constraintValue, dimensions, CompactPolynomialC6(rbfParameter.value),
                                                            {{xDead, yDead, zDead}}, solverRtol, polynomial, preallocation));
    } else if (type == VALUE_AXIAL_GEOMETRIC_MULTISCALE) {

      AxialGeoMultiscaleMapping::MultiscaleType multiscaleTypeValue;
      if (multiscaleType == "spread") {
        multiscaleTypeValue = AxialGeoMultiscaleMapping::SPREAD;
      } else if (multiscaleType == "collect") {
        multiscaleTypeValue = AxialGeoMultiscaleMapping::COLLECT;
      } else {
        PRECICE_ERROR("Unknown geometric multiscale type \"{}\". Known types are \"spread\" and \"collect\".", multiscaleTypeValue);
      }
      configuredMapping.mapping = PtrMapping(
          new AxialGeoMultiscaleMapping(constraintValue, dimensions, multiscaleTypeValue, radius));
      configuredMapping.isRBF = false;
    } else {
      PRECICE_ERROR("Unknown mapping type!");
    }
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
  return mappingType == TYPE_RBF_PUM_DIRECT || mappingType == TYPE_RBF_GLOBAL_DIRECT || mappingType == TYPE_RBF_GLOBAL_ITERATIVE || mappingType == TYPE_RBF_ALIAS;
}

BasisFunction MappingConfiguration::parseBasisFunctions(const std::string &basisFctName) const
{
  BasisFunction basisFunction{};
  if (basisFctName == RBF_TPS)
    basisFunction = BasisFunction::ThinPlateSplines;
  else if (basisFctName == RBF_MULTIQUADRICS)
    basisFunction = BasisFunction::Multiquadrics;
  else if (basisFctName == RBF_INV_MULTIQUADRICS)
    basisFunction = BasisFunction::InverseMultiquadrics;
  else if (basisFctName == RBF_VOLUME_SPLINES)
    basisFunction = BasisFunction::VolumeSplines;
  else if (basisFctName == RBF_GAUSSIAN)
    basisFunction = BasisFunction::Gaussian;
  else if (basisFctName == RBF_CTPS_C2)
    basisFunction = BasisFunction::CompactThinPlateSplinesC2;
  else if (basisFctName == RBF_CPOLYNOMIAL_C0)
    basisFunction = BasisFunction::WendlandC0;
  else if (basisFctName == RBF_CPOLYNOMIAL_C2)
    basisFunction = BasisFunction::WendlandC2;
  else if (basisFctName == RBF_CPOLYNOMIAL_C4)
    basisFunction = BasisFunction::WendlandC4;
  else if (basisFctName == RBF_CPOLYNOMIAL_C6)
    basisFunction = BasisFunction::WendlandC6;
  else
    PRECICE_UNREACHABLE("Unknown basis function \"{}\".", basisFctName);
  return basisFunction;
}
} // namespace precice::mapping
