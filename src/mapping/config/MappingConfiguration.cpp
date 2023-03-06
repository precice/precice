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
using variant_t = std::variant<xml::XMLAttribute<double>, xml::XMLAttribute<std::string>, xml::XMLAttribute<bool>>;
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
  Ginkgo
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

// Specialization for the Ginkgo RBF backend
#ifndef PRECICE_NO_Ginkgo
template <typename RBF>
struct BackendSelector<RBFBackend::Ginkgo, RBF> {
  typedef mapping::RadialBasisFctMapping<RBF> type;
};
#endif

// Variant holding all available RBF classes
using rbf_variant_t = std::variant<CompactPolynomialC0, CompactPolynomialC2, CompactPolynomialC4, CompactPolynomialC6, CompactThinPlateSplinesC2, ThinPlateSplines, VolumeSplines, Multiquadrics, InverseMultiquadrics, Gaussian>;

// The actual instantiation of the mapping class, which is called by the visitor \ref getRBFMapping
template <RBFBackend T, typename RADIAL_BASIS_FUNCTION_T, typename... Args>
PtrMapping instantiateRBFMapping(mapping::Mapping::Constraint &constraint, int dimension, RADIAL_BASIS_FUNCTION_T function,
                                 Args &&... args)
{
  return PtrMapping(new typename BackendSelector<T, RADIAL_BASIS_FUNCTION_T>::type(constraint, dimension, function, std::forward<Args>(args)...));
}

// Constrcuts the RBF function based on the functionType
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

  auto attrSolverRtol = makeXMLAttribute(ATTR_SOLVER_RTOL, 1e-9)
                            .setDocumentation("Solver relative tolerance for convergence");
  auto attrMaxIterations = makeXMLAttribute(ATTR_MAX_ITERATIONS, 1e6)
                               .setDocumentation("Maximum number of iterations of the solver");

  auto attrPreallocation = makeXMLAttribute(ATTR_PREALLOCATION, PREALLOCATION_TREE)
                               .setDocumentation("Sets kind of preallocation for PETSc RBF implementation")
                               .setOptions({PREALLOCATION_ESTIMATE, PREALLOCATION_COMPUTE, PREALLOCATION_OFF, PREALLOCATION_SAVE, PREALLOCATION_TREE});

  auto attrExecutor = makeXMLAttribute(ATTR_EXECUTOR, "reference-executor")
                          .setDocumentation("Specifies the execution backend used by Ginkgo.")
                          .setOptions({"reference-executor", "omp-executor", "cuda-executor", "hip-executor"});

  auto attrDeviceId = makeXMLAttribute(ATTR_DEVICE_ID, static_cast<double>(0))
                          .setDocumentation("Specifies the ID of the GPU that should be used for the Ginkgo GPU backend.");

  auto attrUnifiedMemory = makeXMLAttribute(ATTR_ENABLE_UNIFIED_MEMORY, false)
                               .setDocumentation("If enabled, CUDA Unified Memory will be enabled which allows CUDA to dynamically access RAM.");

  auto attrSolver = makeXMLAttribute(ATTR_SOLVER, "cg-solver")
                        .setDocumentation("Specifies the iterative solver used by Ginkgo.")
                        .setOptions({"cg-solver", "gmres-solver", "mg-solver"});

  auto attrPreconditioner = makeXMLAttribute(ATTR_PRECONDITIONER, "jacobi-preconditioner")
                                .setDocumentation("Specifies the preconditioner used by Ginkgo.")
                                .setOptions({"jacobi-preconditioner", "cholesky-preconditioner", "ilu-preconditioner", "isai-preconditioner", "no-preconditioner"});

  auto attrUsePreconditioner = makeXMLAttribute(ATTR_USE_PRECONDITIONER, true)
                                   .setDocumentation("If enabled, the Ginkgo solver will apply a preconditioner to the linear system")
                                   .setOptions({true, false});

  auto attrJacobiBlockSize = makeXMLAttribute(ATTR_JACOBI_BLOCK_SIZE, static_cast<double>(1)) // TODO: Fix datatype
                                 .setDocumentation("Size of diagonal blocks for Jacobi preconditioner.");

  // Add the relevant attributes to the relevant tags
  addAttributes(projectionTags, {attrFromMesh, attrToMesh, attrDirection, attrConstraint});
  addAttributes(rbfDirectTags, {attrFromMesh, attrToMesh, attrDirection, attrConstraint, attrPolynomial, attrXDead, attrYDead, attrZDead});
  addAttributes(rbfIterativeTags, {attrFromMesh, attrToMesh, attrDirection, attrConstraint, attrPolynomial, attrXDead, attrYDead, attrZDead, attrMaxIterations, attrSolverRtol, attrPreallocation, attrExecutor, attrDeviceId, attrUnifiedMemory, attrSolver, attrUsePreconditioner, attrPreconditioner, attrJacobiBlockSize});
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
  addSubtagsToParents(shapeParameterRBF, rbfAliasTag);

  // For the Gaussian, we need default values as the user can pass a support radius or a shape parameter
  std::list<XMLTag> GaussRBF{
      XMLTag{*this, RBF_GAUSSIAN, once, SUBTAG_BASIS_FUNCTION}.setDocumentation("Gaussian basis function accepting a support radius or a shape parameter.")};
  attrShapeParam.setDefaultValue(std::numeric_limits<double>::quiet_NaN());
  attrSupportRadius.setDefaultValue(std::numeric_limits<double>::quiet_NaN());
  addAttributes(GaussRBF, {attrShapeParam, attrSupportRadius});
  addSubtagsToParents(GaussRBF, rbfIterativeTags);
  addSubtagsToParents(GaussRBF, rbfDirectTags);
  addSubtagsToParents(GaussRBF, rbfAliasTag);

  // tags without an attribute
  std::list<XMLTag> attributelessRBFs{
      XMLTag{*this, RBF_TPS, once, SUBTAG_BASIS_FUNCTION}.setDocumentation("Thin-plate-splines"),
      XMLTag{*this, RBF_VOLUME_SPLINES, once, SUBTAG_BASIS_FUNCTION}.setDocumentation("Volume splines")};

  addSubtagsToParents(attributelessRBFs, rbfIterativeTags);
  addSubtagsToParents(attributelessRBFs, rbfDirectTags);
  addSubtagsToParents(attributelessRBFs, rbfAliasTag);

  // Add all tags to the mapping tag
  parent.addSubtags(projectionTags);
  parent.addSubtags(rbfIterativeTags);
  parent.addSubtags(rbfDirectTags);
  parent.addSubtags(rbfAliasTag);
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
    bool        xDead         = tag.getBooleanAttributeValue(ATTR_X_DEAD, false);
    bool        yDead         = tag.getBooleanAttributeValue(ATTR_Y_DEAD, false);
    bool        zDead         = tag.getBooleanAttributeValue(ATTR_Z_DEAD, false);
    double      solverRtol    = tag.getDoubleAttributeValue(ATTR_SOLVER_RTOL, 1e-9);
    std::string strPolynomial = tag.getStringAttributeValue(ATTR_POLYNOMIAL, POLYNOMIAL_SEPARATE);
    std::string strPrealloc   = tag.getStringAttributeValue(ATTR_PREALLOCATION, PREALLOCATION_TREE);

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

    _rbfConfig = configureRBFMapping(type, context, strPolynomial, strPrealloc, xDead, yDead, zDead, solverRtol);

    _ginkgoParameter.residualNorm        = _rbfConfig.solverRtol;
    _ginkgoParameter.executor            = tag.getStringAttributeValue(ATTR_EXECUTOR, "reference-executor");
    _ginkgoParameter.solver              = tag.getStringAttributeValue(ATTR_SOLVER, "cg-solver");
    _ginkgoParameter.preconditioner      = tag.getStringAttributeValue(ATTR_PRECONDITIONER, "jacobi-preconditioner");
    _ginkgoParameter.usePreconditioner   = tag.getBooleanAttributeValue(ATTR_USE_PRECONDITIONER, true);
    _ginkgoParameter.jacobiBlockSize     = tag.getDoubleAttributeValue(ATTR_JACOBI_BLOCK_SIZE, 4);
    _ginkgoParameter.maxIterations       = tag.getDoubleAttributeValue(ATTR_MAX_ITERATIONS, 1e6);
    _ginkgoParameter.deviceId            = tag.getDoubleAttributeValue(ATTR_DEVICE_ID, 0);
    _ginkgoParameter.enableUnifiedMemory = tag.getBooleanAttributeValue(ATTR_ENABLE_UNIFIED_MEMORY, false);

    checkDuplicates(configuredMapping);
    _mappings.push_back(configuredMapping);
  } else if (tag.getNamespace() == SUBTAG_BASIS_FUNCTION) {

    PRECICE_ASSERT(!_mappings.empty());
    // We can only set one subtag
    PRECICE_CHECK(_mappings.back().mapping == nullptr, "More than one basis-function was defined for the.");

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

    // Instantiate the RBF mapping classes
    if (_rbfConfig.solver == RBFConfiguration::SystemSolver::GlobalDirect) {
      mapping.mapping = getRBFMapping<RBFBackend::Eigen>(basisFunction, constraintValue, mapping.fromMesh->getDimensions(), supportRadius, shapeParameter, _rbfConfig.deadAxis, _rbfConfig.polynomial);
    } else if (_rbfConfig.solver == RBFConfiguration::SystemSolver::GlobalIterative) {
#ifndef PRECICE_NO_PETSC
      // for petsc initialization
      int   argc = 1;
      char *arg  = new char[8];
      strcpy(arg, "precice");
      char **argv = &arg;
      utils::Petsc::initialize(&argc, &argv, utils::Parallel::current()->comm);
      delete[] arg;

      mapping.mapping = getRBFMapping<RBFBackend::PETSc>(basisFunction, constraintValue, mapping.fromMesh->getDimensions(), supportRadius, shapeParameter, _rbfConfig.deadAxis, _rbfConfig.solverRtol, _rbfConfig.polynomial, _rbfConfig.preallocation);

#elif !defined(PRECICE_NO_GINKGO)
      mapping.mapping = getRBFMapping<RBFBackend::Ginkgo>(basisFunction, constraintValue, mapping.fromMesh->getDimensions(), supportRadius, shapeParameter, _rbfConfig.deadAxis, _rbfConfig.polynomial, false, _ginkgoParameter);

#else
      PRECICE_CHECK(false, "The global-iterative RBF solver requires a preCICE build with PETSc enabled.");
#endif
    } else {
      PRECICE_UNREACHABLE("Unknown RBF solver.");
    }
  }
}

MappingConfiguration::RBFConfiguration MappingConfiguration::configureRBFMapping(const std::string &              type,
                                                                                 const xml::ConfigurationContext &context,
                                                                                 const std::string &              polynomial,
                                                                                 const std::string &              preallocation,
                                                                                 bool xDead, bool yDead, bool zDead,
                                                                                 double solverRtol) const
{
  RBFConfiguration rbfConfig;

  if (type == TYPE_RBF_GLOBAL_ITERATIVE)
    rbfConfig.solver = RBFConfiguration::SystemSolver::GlobalIterative;
  else if (type == TYPE_RBF_GLOBAL_DIRECT)
    rbfConfig.solver = RBFConfiguration::SystemSolver::GlobalDirect;
  else {
    // Rather simple auto-selection (for now)
    // The default is the Eigen backend, as it is always available. Only in certain situations, we will decide for the PETSc backend
    rbfConfig.solver = RBFConfiguration::SystemSolver::GlobalDirect;

#ifndef PRECICE_NO_PETSC
    // Running in serial, the Eigen backend will most likely be the fastest and most accurate variant
    // We decide for the PETSc variant only if we have a lot of interface vertices and a lot of ranks
    // A more sophisticated criterion here could take the globalNumberOfVertices into account
    // (the mesh pointer is stored in the configuredMapping anyway), but this quantity is not yet computed
    // during the configuration time.
    if (context.size > 16) {
      rbfConfig.solver = RBFConfiguration::SystemSolver::GlobalIterative;
    }
#endif
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
    PRECICE_UNREACHABLE("Unknwon preallocation configuration");

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
  return mappingType == TYPE_RBF_GLOBAL_DIRECT || mappingType == TYPE_RBF_GLOBAL_ITERATIVE || mappingType == TYPE_RBF_ALIAS;
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
