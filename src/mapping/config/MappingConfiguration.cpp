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
#include "mapping/GinkgoRadialBasisFctSolver.hpp"
#include "mapping/LinearCellInterpolationMapping.hpp"
#include "mapping/Mapping.hpp"
#include "mapping/NearestNeighborGradientMapping.hpp"
#include "mapping/NearestNeighborMapping.hpp"
#include "mapping/NearestProjectionMapping.hpp"
#include "mapping/PartitionOfUnityMapping.hpp"
#include "mapping/PetRadialBasisFctMapping.hpp"
#include "mapping/RadialBasisFctMapping.hpp"
#include "mapping/RadialBasisFctSolver.hpp"
#include "mapping/RadialGeoMultiscaleMapping.hpp"
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
  for (auto &p : parents) {
    p.addSubtags(subtags);
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
  Ginkgo,
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
  typedef mapping::RadialBasisFctMapping<RadialBasisFctSolver<RBF>> type;
};

// Specialization for the PETSc RBF backend
#ifndef PRECICE_NO_PETSC
template <typename RBF>
struct BackendSelector<RBFBackend::PETSc, RBF> {
  typedef mapping::PetRadialBasisFctMapping<RBF> type;
};
#endif

// Specialization for the Ginkgo RBF backend
#ifndef PRECICE_NO_GINKGO
template <typename RBF>
struct BackendSelector<RBFBackend::Ginkgo, RBF> {
  typedef mapping::RadialBasisFctMapping<GinkgoRadialBasisFctSolver<RBF>, MappingConfiguration::GinkgoParameter> type;
};
#endif
// Specialization for the RBF PUM backend
template <typename RBF>
struct BackendSelector<RBFBackend::PUM, RBF> {
  typedef mapping::PartitionOfUnityMapping<RBF> type;
};

// Variant holding all available RBF classes
using rbf_variant_t = std::variant<CompactPolynomialC0, CompactPolynomialC2, CompactPolynomialC4, CompactPolynomialC6, CompactPolynomialC8, CompactThinPlateSplinesC2, ThinPlateSplines, VolumeSplines, Multiquadrics, InverseMultiquadrics, Gaussian>;

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
  case BasisFunction::WendlandC8: {
    return mapping::CompactPolynomialC8(supportRadius);
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
  std::list<XMLTag> geoMultiscaleTags{
      XMLTag{*this, TYPE_AXIAL_GEOMETRIC_MULTISCALE, occ, TAG}.setDocumentation("Axial geometric multiscale mapping between one 1D and multiple 3D vertices."),
      XMLTag{*this, TYPE_RADIAL_GEOMETRIC_MULTISCALE, occ, TAG}.setDocumentation("Radial geometric multiscale mapping between multiple 1D and multiple 3D vertices, distributed along a principle axis.")};

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
  // TODO: Discuss whether we wanto to introduce this attribute
  // auto attrMaxIterations = makeXMLAttribute(ATTR_MAX_ITERATIONS, 1e6)
  //                              .setDocumentation("Maximum number of iterations of the solver");

  auto verticesPerCluster = XMLAttribute<int>(ATTR_VERTICES_PER_CLUSTER, 50)
                                .setDocumentation("Average number of vertices per cluster (partition) applied in the rbf partition of unity method.");
  auto relativeOverlap = makeXMLAttribute(ATTR_RELATIVE_OVERLAP, 0.15)
                             .setDocumentation("Value between 0 and 1 indicating the relative overlap between clusters. A value of 0.15 is usually a good trade-off between accuracy and efficiency.");
  auto projectToInput = XMLAttribute<bool>(ATTR_PROJECT_TO_INPUT, true)
                            .setDocumentation("If enabled, places the cluster centers at the closest vertex of the input mesh. Should be enabled in case of non-uniform point distributions such as for shell structures.");

  auto attrGeoMultiscaleType = XMLAttribute<std::string>(ATTR_GEOMETRIC_MULTISCALE_TYPE)
                                   .setDocumentation("Type of geometric multiscale mapping. Either 'spread' or 'collect'.")
                                   .setOptions({GEOMETRIC_MULTISCALE_TYPE_SPREAD, GEOMETRIC_MULTISCALE_TYPE_COLLECT});
  auto attrGeoMultiscaleAxis = XMLAttribute<std::string>(ATTR_GEOMETRIC_MULTISCALE_AXIS)
                                   .setDocumentation("Principle axis along which geometric multiscale mapping is performed.")
                                   .setOptions({GEOMETRIC_MULTISCALE_AXIS_X, GEOMETRIC_MULTISCALE_AXIS_Y, GEOMETRIC_MULTISCALE_AXIS_Z});
  auto attrGeoMultiscaleRadius = XMLAttribute<double>(ATTR_GEOMETRIC_MULTISCALE_RADIUS)
                                     .setDocumentation("Radius of the circular interface between the 1D and 3D participant.");

  // Add the relevant attributes to the relevant tags
  addAttributes(projectionTags, {attrFromMesh, attrToMesh, attrDirection, attrConstraint});
  addAttributes(rbfDirectTags, {attrFromMesh, attrToMesh, attrDirection, attrConstraint, attrPolynomial, attrXDead, attrYDead, attrZDead});
  addAttributes(rbfIterativeTags, {attrFromMesh, attrToMesh, attrDirection, attrConstraint, attrPolynomial, attrXDead, attrYDead, attrZDead, attrSolverRtol});
  addAttributes(pumDirectTags, {attrFromMesh, attrToMesh, attrDirection, attrConstraint, attrPumPolynomial, verticesPerCluster, relativeOverlap, projectToInput});
  addAttributes(rbfAliasTag, {attrFromMesh, attrToMesh, attrDirection, attrConstraint, attrXDead, attrYDead, attrZDead});
  addAttributes(geoMultiscaleTags, {attrFromMesh, attrToMesh, attrDirection, attrConstraint, attrGeoMultiscaleType, attrGeoMultiscaleAxis, attrGeoMultiscaleRadius});

  // Now we take care of the subtag executor. We repeat some of the subtags in order to add individual documentation
  XMLTag::Occurrence once = XMLTag::OCCUR_NOT_OR_ONCE;
  // TODO, make type an int
  auto attrDeviceId = makeXMLAttribute(ATTR_DEVICE_ID, static_cast<int>(0))
                          .setDocumentation("Specifies the ID of the GPU that should be used for the Ginkgo GPU backend.");
  auto attrNThreads = makeXMLAttribute(ATTR_N_THREADS, static_cast<int>(0))
                          .setDocumentation("Specifies the number of threads for the OpenMP executor that should be used for the Ginkgo OpenMP backend. If a value of \"0\" is set, preCICE doesn't set the number of threads and the default behavior of OpenMP applies.");

  // First, we have the executors for the direct solvers
  {
    std::list<XMLTag> cpuExecutor{
        XMLTag{*this, EXECUTOR_CPU, once, SUBTAG_EXECUTOR}.setDocumentation("The default executor, which uses a single-core CPU with a gather-scatter parallelism.")};
    std::list<XMLTag> deviceExecutors{
        XMLTag{*this, EXECUTOR_CUDA, once, SUBTAG_EXECUTOR}.setDocumentation("Cuda (Nvidia) executor, which uses cuSolver/Ginkgo and a direct QR decomposition with a gather-scatter parallelism."),
        XMLTag{*this, EXECUTOR_HIP, once, SUBTAG_EXECUTOR}.setDocumentation("Hip (AMD/Nvidia) executor, which uses hipSolver/Ginkgo and a direct QR decomposition with a gather-scatter parallelism.")};

    addAttributes(deviceExecutors, {attrDeviceId});
    addSubtagsToParents(cpuExecutor, rbfDirectTags);
    addSubtagsToParents(deviceExecutors, rbfDirectTags);
  }
  // Second, the executors for the iterative solver
  {
    std::list<XMLTag> cpuExecutor{
        XMLTag{*this, EXECUTOR_CPU, once, SUBTAG_EXECUTOR}.setDocumentation("The default executor relying on PETSc, which uses CPUs and distributed memory parallelism via MPI.")};
    std::list<XMLTag> deviceExecutors{
        XMLTag{*this, EXECUTOR_CUDA, once, SUBTAG_EXECUTOR}.setDocumentation("Cuda (Nvidia) executor, which uses Ginkgo with a gather-scatter parallelism."),
        XMLTag{*this, EXECUTOR_HIP, once, SUBTAG_EXECUTOR}.setDocumentation("Hip (AMD/Nvidia) executor, which uses hipSolver with a gather-scatter parallelism.")};
    std::list<XMLTag> ompExecutor{
        XMLTag{*this, EXECUTOR_OMP, once, SUBTAG_EXECUTOR}.setDocumentation("OpenMP executor, which uses Ginkgo with a gather-scatter parallelism.")};

    addAttributes(deviceExecutors, {attrDeviceId});
    addAttributes(ompExecutor, {attrNThreads});
    addSubtagsToParents(cpuExecutor, rbfIterativeTags);
    addSubtagsToParents(deviceExecutors, rbfIterativeTags);
    addSubtagsToParents(ompExecutor, rbfIterativeTags);
  }
  {
    std::list<XMLTag> cpuExecutor{
        XMLTag{*this, EXECUTOR_CPU, once, SUBTAG_EXECUTOR}.setDocumentation("The default (and currently only) executor using a CPU and a distributed memory parallelism via MPI.")};
    addSubtagsToParents(cpuExecutor, pumDirectTags);
  }
  // The alias tag doesn't receive the subtag at all

  // Now we take care of the subtag basis function
  // First, we have the tags using a support radius
  std::list<XMLTag> supportRadiusRBF{
      XMLTag{*this, RBF_CPOLYNOMIAL_C0, once, SUBTAG_BASIS_FUNCTION}.setDocumentation("Wendland C0 function"),
      XMLTag{*this, RBF_CPOLYNOMIAL_C2, once, SUBTAG_BASIS_FUNCTION}.setDocumentation("Wendland C2 function"),
      XMLTag{*this, RBF_CPOLYNOMIAL_C4, once, SUBTAG_BASIS_FUNCTION}.setDocumentation("Wendland C4 function"),
      XMLTag{*this, RBF_CPOLYNOMIAL_C6, once, SUBTAG_BASIS_FUNCTION}.setDocumentation("Wendland C6 function"),
      XMLTag{*this, RBF_CPOLYNOMIAL_C8, once, SUBTAG_BASIS_FUNCTION}.setDocumentation("Wendland C8 function"),
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
  parent.addSubtags(geoMultiscaleTags);
}

void MappingConfiguration::setExperimental(
    bool experimental)
{
  _experimental = experimental;
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
    // We set here default values, but their actual value doesn't really matter.
    // It's just for the mapping methods, which do not use these attributes at all.
    bool        xDead         = tag.getBooleanAttributeValue(ATTR_X_DEAD, false);
    bool        yDead         = tag.getBooleanAttributeValue(ATTR_Y_DEAD, false);
    bool        zDead         = tag.getBooleanAttributeValue(ATTR_Z_DEAD, false);
    double      solverRtol    = tag.getDoubleAttributeValue(ATTR_SOLVER_RTOL, 1e-9);
    std::string strPolynomial = tag.getStringAttributeValue(ATTR_POLYNOMIAL, POLYNOMIAL_SEPARATE);

    // geometric multiscale related tags
    std::string geoMultiscaleType = tag.getStringAttributeValue(ATTR_GEOMETRIC_MULTISCALE_TYPE, "");
    std::string geoMultiscaleAxis = tag.getStringAttributeValue(ATTR_GEOMETRIC_MULTISCALE_AXIS, "");
    double      multiscaleRadius  = tag.getDoubleAttributeValue(ATTR_GEOMETRIC_MULTISCALE_RADIUS, 1.0);

    if (type == TYPE_AXIAL_GEOMETRIC_MULTISCALE || type == TYPE_RADIAL_GEOMETRIC_MULTISCALE) {
      PRECICE_CHECK(_experimental, "Axial geometric multiscale is experimental and the configuration can change between minor releases. Set experimental=\"on\" in the precice-configuration tag.");
    }

    if (type == TYPE_AXIAL_GEOMETRIC_MULTISCALE && context.size > 1) {
      throw std::runtime_error{"Axial geometric multiscale mapping is not available for parallel participants."};
    }

    if (type == TYPE_RADIAL_GEOMETRIC_MULTISCALE && context.size > 1) {
      throw std::runtime_error{"Radial geometric multiscale mapping is not available for parallel participants."};
    }

    // pum related tags
    int    verticesPerCluster = tag.getIntAttributeValue(ATTR_VERTICES_PER_CLUSTER, 100);
    double relativeOverlap    = tag.getDoubleAttributeValue(ATTR_RELATIVE_OVERLAP, 0.3);
    bool   projectToInput     = tag.getBooleanAttributeValue(ATTR_PROJECT_TO_INPUT, true);

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

    ConfiguredMapping configuredMapping = createMapping(dir, type, fromMesh, toMesh, geoMultiscaleType, geoMultiscaleAxis, multiscaleRadius);

    _rbfConfig = configureRBFMapping(type, strPolynomial, xDead, yDead, zDead, solverRtol, verticesPerCluster, relativeOverlap, projectToInput);

    checkDuplicates(configuredMapping);
    _mappings.push_back(configuredMapping);
  } else if (tag.getNamespace() == SUBTAG_BASIS_FUNCTION) {

    PRECICE_ASSERT(!_mappings.empty());
    PRECICE_CHECK(_mappings.back().requiresBasisFunction == true, "A basis-function was defined for the mapping "
                                                                  "from mesh \"{}\" to mesh \"{}\", but no basis-function is required for this mapping type. "
                                                                  "Please remove the basis-function tag or configure an rbf mapping, which requires a basis-function.",
                  _mappings.back().fromMesh->getName(), _mappings.back().toMesh->getName());
    // We can only set one subtag
    PRECICE_CHECK(_rbfConfig.basisFunctionDefined == false, "More than one basis-function was defined for the mapping "
                                                            "from mesh \"{}\" to mesh \"{}\".",
                  _mappings.back().fromMesh->getName(), _mappings.back().toMesh->getName());

    std::string basisFctName   = tag.getName();
    double      supportRadius  = tag.getDoubleAttributeValue(ATTR_SUPPORT_RADIUS, 0.);
    double      shapeParameter = tag.getDoubleAttributeValue(ATTR_SHAPE_PARAM, 0.);

    _rbfConfig.basisFunction        = parseBasisFunctions(basisFctName);
    _rbfConfig.basisFunctionDefined = true;
    // The Gaussian RBF is always treated as a shape-parameter RBF. Hence, we have to convert the support radius, if necessary
    if (_rbfConfig.basisFunction == BasisFunction::Gaussian) {
      const bool exactlyOneSet = (std::isfinite(supportRadius) && !std::isfinite(shapeParameter)) ||
                                 (std::isfinite(shapeParameter) && !std::isfinite(supportRadius));
      PRECICE_CHECK(exactlyOneSet, "The specified parameters for the Gaussian RBF mapping are invalid. Please specify either a \"shape-parameter\" or a \"support-radius\".");

      if (std::isfinite(supportRadius) && !std::isfinite(shapeParameter)) {
        shapeParameter = std::sqrt(-std::log(Gaussian::cutoffThreshold)) / supportRadius;
      }
    }

    _rbfConfig.supportRadius  = supportRadius;
    _rbfConfig.shapeParameter = shapeParameter;
  } else if (tag.getNamespace() == SUBTAG_EXECUTOR) {
    _executorConfig = std::make_unique<ExecutorConfiguration>();

    if (tag.getName() == EXECUTOR_CPU) {
      _executorConfig->executor = ExecutorConfiguration::Executor::CPU;
    } else if (tag.getName() == EXECUTOR_CUDA) {
      _executorConfig->executor = ExecutorConfiguration::Executor::CUDA;
    } else if (tag.getName() == EXECUTOR_HIP) {
      _executorConfig->executor = ExecutorConfiguration::Executor::HIP;
    } else if (tag.getName() == EXECUTOR_OMP) {
      _executorConfig->executor = ExecutorConfiguration::Executor::OpenMP;
    }

    _executorConfig->deviceId = tag.getIntAttributeValue(ATTR_DEVICE_ID, -1);
    _executorConfig->nThreads = tag.getIntAttributeValue(ATTR_N_THREADS, 0);
  }
}

MappingConfiguration::RBFConfiguration MappingConfiguration::configureRBFMapping(const std::string &type,
                                                                                 const std::string &polynomial,
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

  rbfConfig.deadAxis   = {{xDead, yDead, zDead}};
  rbfConfig.solverRtol = solverRtol;

  rbfConfig.verticesPerCluster = verticesPerCluster;
  rbfConfig.relativeOverlap    = relativeOverlap;
  rbfConfig.projectToInput     = projectToInput;

  return rbfConfig;
}

MappingConfiguration::ConfiguredMapping MappingConfiguration::createMapping(
    const std::string &direction,
    const std::string &type,
    const std::string &fromMeshName,
    const std::string &toMeshName,
    const std::string &geoMultiscaleType,
    const std::string &geoMultiscaleAxis,
    const double &     multiscaleRadius) const
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

  // Check for compatible mesh dimensions
  PRECICE_CHECK(fromMesh->getDimensions() == toMesh->getDimensions(),
                "Mapping between meshes of different dimensions is not allowed yet. "
                "Please set the same dimensions attribute to meshes \"{}\" and \"{}\", "
                "or choose different meshes.",
                fromMeshName, toMeshName);

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

  } else if (type == TYPE_AXIAL_GEOMETRIC_MULTISCALE) {

    // Axial geometric multiscale is not applicable with the conservative constraint
    PRECICE_CHECK(constraintValue != Mapping::CONSERVATIVE,
                  "Axial geometric multiscale mapping is not implemented for the \"conservative\" constraint. "
                  "Please select constraint=\" consistent\" or a different mapping method.");

    // Convert strings into enums
    AxialGeoMultiscaleMapping::MultiscaleAxis multiscaleAxis;
    if (geoMultiscaleAxis == "x") {
      multiscaleAxis = AxialGeoMultiscaleMapping::MultiscaleAxis::X;
    } else if (geoMultiscaleAxis == "y") {
      multiscaleAxis = AxialGeoMultiscaleMapping::MultiscaleAxis::Y;
    } else if (geoMultiscaleAxis == "z") {
      multiscaleAxis = AxialGeoMultiscaleMapping::MultiscaleAxis::Z;
    } else {
      PRECICE_UNREACHABLE("Unknown geometric multiscale axis \"{}\".", geoMultiscaleAxis);
    }

    AxialGeoMultiscaleMapping::MultiscaleType multiscaleType;
    if (geoMultiscaleType == "spread") {
      multiscaleType = AxialGeoMultiscaleMapping::MultiscaleType::SPREAD;
    } else if (geoMultiscaleType == "collect") {
      multiscaleType = AxialGeoMultiscaleMapping::MultiscaleType::COLLECT;
    } else {
      PRECICE_UNREACHABLE("Unknown geometric multiscale type \"{}\".", geoMultiscaleType);
    }

    configuredMapping.mapping = PtrMapping(new AxialGeoMultiscaleMapping(constraintValue, fromMesh->getDimensions(), multiscaleType, multiscaleAxis, multiscaleRadius));

  } else if (type == TYPE_RADIAL_GEOMETRIC_MULTISCALE) {

    // Radial geometric multiscale is not applicable with the conservative constraint
    PRECICE_CHECK(constraintValue != Mapping::CONSERVATIVE,
                  "Radial geometric multiscale mapping is not implemented for the \"conservative\" constraint. "
                  "Please select constraint=\" consistent\" or a different mapping method.");

    // Convert strings into enums
    RadialGeoMultiscaleMapping::MultiscaleAxis multiscaleAxis;
    if (geoMultiscaleAxis == "x") {
      multiscaleAxis = RadialGeoMultiscaleMapping::MultiscaleAxis::X;
    } else if (geoMultiscaleAxis == "y") {
      multiscaleAxis = RadialGeoMultiscaleMapping::MultiscaleAxis::Y;
    } else if (geoMultiscaleAxis == "z") {
      multiscaleAxis = RadialGeoMultiscaleMapping::MultiscaleAxis::Z;
    } else {
      PRECICE_UNREACHABLE("Unknown geometric multiscale axis \"{}\".", geoMultiscaleAxis);
    }

    RadialGeoMultiscaleMapping::MultiscaleType multiscaleType;
    if (geoMultiscaleType == "spread") {
      multiscaleType = RadialGeoMultiscaleMapping::MultiscaleType::SPREAD;
    } else if (geoMultiscaleType == "collect") {
      multiscaleType = RadialGeoMultiscaleMapping::MultiscaleType::COLLECT;
    } else {
      PRECICE_UNREACHABLE("Unknown geometric multiscale type \"{}\".", geoMultiscaleType);
    }

    configuredMapping.mapping = PtrMapping(new RadialGeoMultiscaleMapping(constraintValue, fromMesh->getDimensions(), multiscaleType, multiscaleAxis));

  } else {
    // We need knowledge about the basis function in order to instantiate the rbf related mapping
    PRECICE_ASSERT(requiresBasisFunction(type));
    configuredMapping.mapping                = nullptr;
    configuredMapping.configuredWithAliasTag = type == TYPE_RBF_ALIAS;
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
  if (tag.getNamespace() == TAG) {
    if (requiresBasisFunction(tag.getName())) {
      PRECICE_CHECK(_rbfConfig.basisFunctionDefined, "No basis-function was defined for the \"{}\" mapping from mesh \"{}\" to mesh \"{}\".", tag.getName(), _mappings.back().fromMesh->getName(), _mappings.back().toMesh->getName());
      if (!_executorConfig) {
        _executorConfig = std::make_unique<ExecutorConfiguration>();
      }
      finishRBFConfiguration();
      _executorConfig.reset();
    }
    PRECICE_ASSERT(_mappings.back().mapping != nullptr);
  }
}

void MappingConfiguration::finishRBFConfiguration()
{
  PRECICE_ASSERT(_executorConfig);
  ConfiguredMapping &mapping = _mappings.back();
  // Instantiate the RBF mapping classes
  // We first categorize according to the executor
  // 1. the CPU executor
  if (_executorConfig->executor == ExecutorConfiguration::Executor::CPU) {
    if (_rbfConfig.solver == RBFConfiguration::SystemSolver::GlobalDirect) {
      mapping.mapping = getRBFMapping<RBFBackend::Eigen>(_rbfConfig.basisFunction, constraintValue, mapping.fromMesh->getDimensions(), _rbfConfig.supportRadius, _rbfConfig.shapeParameter, _rbfConfig.deadAxis, _rbfConfig.polynomial);
    } else if (_rbfConfig.solver == RBFConfiguration::SystemSolver::GlobalIterative) {
#ifndef PRECICE_NO_PETSC
      // for petsc initialization
      utils::Petsc::initialize(utils::Parallel::current()->comm);

      mapping.mapping = getRBFMapping<RBFBackend::PETSc>(_rbfConfig.basisFunction, constraintValue, mapping.fromMesh->getDimensions(), _rbfConfig.supportRadius, _rbfConfig.shapeParameter, _rbfConfig.deadAxis, _rbfConfig.solverRtol, _rbfConfig.polynomial);
#else
      PRECICE_CHECK(false, "The global-iterative RBF solver on a CPU requires a preCICE build with PETSc enabled.");
#endif
    } else if (_rbfConfig.solver == RBFConfiguration::SystemSolver::PUMDirect) {
      mapping.mapping = getRBFMapping<RBFBackend::PUM>(_rbfConfig.basisFunction, constraintValue, mapping.fromMesh->getDimensions(), _rbfConfig.supportRadius, _rbfConfig.shapeParameter, _rbfConfig.polynomial, _rbfConfig.verticesPerCluster, _rbfConfig.relativeOverlap, _rbfConfig.projectToInput);
    } else {
      PRECICE_UNREACHABLE("Unknown RBF solver.");
    }
    // 2. any other executor is configured via Ginkgo
  } else {
#ifndef PRECICE_NO_GINKGO
    _ginkgoParameter                   = GinkgoParameter();
    _ginkgoParameter.usePreconditioner = false;
    _ginkgoParameter.deviceId          = _executorConfig->deviceId;
    if (_executorConfig->executor == ExecutorConfiguration::Executor::CUDA) {
      _ginkgoParameter.executor = "cuda-executor";
#ifndef PRECICE_WITH_CUDA
      PRECICE_CHECK(false, "The cuda-executor (configured for the mapping from mesh {} to mesh {}) requires a Ginkgo and preCICE build with Cuda enabled.", mapping.fromMesh->getName(), mapping.toMesh->getName());
#endif
    } else if (_executorConfig->executor == ExecutorConfiguration::Executor::HIP) {
      _ginkgoParameter.executor = "hip-executor";
#ifndef PRECICE_WITH_HIP
      PRECICE_CHECK(false, "The hip-executor (configured for the mapping from mesh {} to mesh {}) requires a Ginkgo and preCICE build with HIP enabled.", mapping.fromMesh->getName(), mapping.toMesh->getName());
#endif
    } else if (_executorConfig->executor == ExecutorConfiguration::Executor::OpenMP) {
      _ginkgoParameter.executor = "omp-executor";
      _ginkgoParameter.nThreads = _executorConfig->nThreads;
#ifndef PRECICE_WITH_OMP
      PRECICE_CHECK(false, "The omp-executor (configured for the mapping from mesh {} to mesh {}) requires a Ginkgo and preCICE build with OpenMP enabled.", mapping.fromMesh->getName(), mapping.toMesh->getName());
#endif
    }
    if (_rbfConfig.solver == RBFConfiguration::SystemSolver::GlobalDirect) {
      _ginkgoParameter.solver = "qr-solver";
    } else if (_rbfConfig.solver == RBFConfiguration::SystemSolver::GlobalIterative) {
      _ginkgoParameter.solver       = "cg-solver";
      _ginkgoParameter.residualNorm = _rbfConfig.solverRtol;
    } else {
      PRECICE_UNREACHABLE("Unknown solver type.");
    }
    mapping.mapping = getRBFMapping<RBFBackend::Ginkgo>(_rbfConfig.basisFunction, constraintValue, mapping.fromMesh->getDimensions(), _rbfConfig.supportRadius, _rbfConfig.shapeParameter, _rbfConfig.deadAxis, _rbfConfig.polynomial, _ginkgoParameter);
#else
    PRECICE_CHECK(false, "The selected executor for the mapping from mesh {} to mesh {} requires a preCICE build with Ginkgo enabled.", mapping.fromMesh->getName(), mapping.toMesh->getName());
#endif
  }
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
  else if (basisFctName == RBF_CPOLYNOMIAL_C8)
    basisFunction = BasisFunction::WendlandC8;
  else
    PRECICE_UNREACHABLE("Unknown basis function \"{}\".", basisFctName);
  return basisFunction;
}
} // namespace precice::mapping
