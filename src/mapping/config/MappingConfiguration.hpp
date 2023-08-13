#pragma once

#include <string>
#include <vector>
#include "logging/Logger.hpp"
#include "mapping/Mapping.hpp"
#include "mapping/SharedPointer.hpp"
#include "mapping/config/MappingConfigurationTypes.hpp"
#include "mesh/SharedPointer.hpp"
#include "xml/XMLTag.hpp"

namespace precice::mapping {

/// Performs XML configuration and holds configured mappings.
class MappingConfiguration : public xml::XMLTag::Listener {
public:
  /// Constants defining the direction of a mapping.
  enum Direction {
    WRITE,
    READ
  };

  /// Configuration data for one mapping.
  struct ConfiguredMapping {
    /// Mapping object.
    PtrMapping mapping;
    /// Remote mesh to map from
    mesh::PtrMesh fromMesh;
    /// Remote mesh to map to
    mesh::PtrMesh toMesh;
    /// Direction of mapping (important to set input and output mesh).
    Direction direction;
    /// true for RBF mapping
    bool requiresBasisFunction;
    /// used the automatic rbf alias tag in order to set the mapping
    bool configuredWithAliasTag = false;
  };

  struct GinkgoParameter {
    std::string  executor            = "reference-executor";
    std::string  solver              = "cg-solver";
    std::string  preconditioner      = "jacobi-preconditioner";
    double       residualNorm        = 1e-8;
    std::size_t  maxIterations       = 1e6;
    bool         usePreconditioner   = false;
    unsigned int jacobiBlockSize     = 4;
    unsigned int deviceId            = 0;
    unsigned int nThreads            = 0;
    bool         enableUnifiedMemory = false;
  };

  MappingConfiguration(
      xml::XMLTag &              parent,
      mesh::PtrMeshConfiguration meshConfiguration);

  /**
   * @brief Callback function required for use of automatic configuration.
   *
   * @return True, if successful.
   */
  virtual void xmlTagCallback(
      const xml::ConfigurationContext &context,
      xml::XMLTag &                    callingTag);

  /**
   * @brief Callback function required for use of automatic configuration.
   *
   * @return True, if successful.
   */
  virtual void xmlEndTagCallback(
      const xml::ConfigurationContext &context,
      xml::XMLTag &                    callingTag);

  /// Returns all configured mappings.
  const std::vector<ConfiguredMapping> &mappings();

  // Only relevant for RBF related mappings
  // Being public here is only required for testing purposes
  struct RBFConfiguration {

    enum struct SystemSolver {
      GlobalDirect,
      GlobalIterative,
      PUMDirect
    };
    SystemSolver        solver{};
    std::array<bool, 3> deadAxis{};
    Polynomial          polynomial{};
    double              solverRtol{};
    int                 verticesPerCluster{};
    double              relativeOverlap{};
    bool                projectToInput{};
    BasisFunction       basisFunction{};
    double              supportRadius{};
    double              shapeParameter{};
    bool                basisFunctionDefined = false;
  };

  struct GeoMultiscaleConfiguration {

    double      multiscaleRadius{};
    std::string multiscaleType{};
    int         multiscaleAxis{};
  };

  /// Returns the RBF configuration, which was configured at the latest.
  /// Only required for the configuration test
  const RBFConfiguration &rbfConfig() const
  {
    return _rbfConfig;
  }

  void resetMappings()
  {
    _mappings.clear();
  }

private:
  mutable logging::Logger _log{"config:MappingConfiguration"};

  const std::string TAG = "mapping";

  // First, declare common attributes and associated options
  const std::string ATTR_TYPE                        = "type";
  const std::string TYPE_NEAREST_NEIGHBOR            = "nearest-neighbor";
  const std::string TYPE_NEAREST_NEIGHBOR_GRADIENT   = "nearest-neighbor-gradient";
  const std::string TYPE_NEAREST_PROJECTION          = "nearest-projection";
  const std::string TYPE_LINEAR_CELL_INTERPOLATION   = "linear-cell-interpolation";
  const std::string TYPE_RBF_GLOBAL_DIRECT           = "rbf-global-direct";
  const std::string TYPE_RBF_GLOBAL_ITERATIVE        = "rbf-global-iterative";
  const std::string TYPE_RBF_PUM_DIRECT              = "rbf-pum-direct";
  const std::string TYPE_RBF_ALIAS                   = "rbf";
  const std::string TYPE_AXIAL_GEOMETRIC_MULTISCALE  = "axial-geometric-multiscale";
  const std::string TYPE_RADIAL_GEOMETRIC_MULTISCALE = "radial-geometric-multiscale";

  const std::string ATTR_DIRECTION  = "direction";
  const std::string DIRECTION_WRITE = "write";
  const std::string DIRECTION_READ  = "read";

  const std::string ATTR_FROM = "from";
  const std::string ATTR_TO   = "to";

  const std::string ATTR_CONSTRAINT                      = "constraint";
  const std::string CONSTRAINT_CONSISTENT                = "consistent";
  const std::string CONSTRAINT_CONSERVATIVE              = "conservative";
  const std::string CONSTRAINT_SCALED_CONSISTENT_SURFACE = "scaled-consistent-surface";
  const std::string CONSTRAINT_SCALED_CONSISTENT_VOLUME  = "scaled-consistent-volume";

  // RBF specific options
  const std::string ATTR_X_DEAD = "x-dead";
  const std::string ATTR_Y_DEAD = "y-dead";
  const std::string ATTR_Z_DEAD = "z-dead";

  const std::string ATTR_POLYNOMIAL     = "polynomial";
  const std::string POLYNOMIAL_SEPARATE = "separate";
  const std::string POLYNOMIAL_ON       = "on";
  const std::string POLYNOMIAL_OFF      = "off";

  // For iterative RBFs
  const std::string ATTR_SOLVER_RTOL = "solver-rtol";

  // For the future
  // const std::string ATTR_PARALLELISM           = "parallelism";
  // const std::string PARALLELISM_GATHER_SCATTER = "gather-scatter";
  // const std::string PARALLELISM                = "distributed";

  // For PUM
  const std::string ATTR_VERTICES_PER_CLUSTER = "vertices-per-cluster";
  const std::string ATTR_RELATIVE_OVERLAP     = "relative-overlap";
  const std::string ATTR_PROJECT_TO_INPUT     = "project-to-input";

  // We declare the basis function as subtag
  const std::string SUBTAG_BASIS_FUNCTION = "basis-function";
  const std::string RBF_TPS               = "thin-plate-splines";
  const std::string RBF_MULTIQUADRICS     = "multiquadrics";
  const std::string RBF_INV_MULTIQUADRICS = "inverse-multiquadrics";
  const std::string RBF_VOLUME_SPLINES    = "volume-splines";
  const std::string RBF_GAUSSIAN          = "gaussian";
  const std::string RBF_CTPS_C2           = "compact-tps-c2";
  const std::string RBF_CPOLYNOMIAL_C0    = "compact-polynomial-c0";
  const std::string RBF_CPOLYNOMIAL_C2    = "compact-polynomial-c2";
  const std::string RBF_CPOLYNOMIAL_C4    = "compact-polynomial-c4";
  const std::string RBF_CPOLYNOMIAL_C6    = "compact-polynomial-c6";

  // Attributes for the subtag
  const std::string ATTR_SHAPE_PARAM    = "shape-parameter";
  const std::string ATTR_SUPPORT_RADIUS = "support-radius";

  // Attributes for geometric multiscale
  const std::string ATTR_GEOMETRIC_MULTISCALE_TYPE   = "multiscale-type";
  const std::string ATTR_GEOMETRIC_MULTISCALE_AXIS   = "multiscale-axis";
  const std::string ATTR_GEOMETRIC_MULTISCALE_RADIUS = "multiscale-radius";
  // For iterative RBFs using Ginkgo
  const std::string SUBTAG_EXECUTOR = "executor";
  const std::string EXECUTOR_CPU    = "cpu";
  const std::string EXECUTOR_CUDA   = "cuda";
  const std::string EXECUTOR_HIP    = "hip";
  const std::string EXECUTOR_OMP    = "openmp";

  const std::string ATTR_DEVICE_ID = "gpu-device-id";
  const std::string ATTR_N_THREADS = "n-threads";
  // const std::string ATTR_ENABLE_UNIFIED_MEMORY = "enable-unified-memory";
  // const std::string ATTR_SOLVER                = "solver";
  // const std::string ATTR_USE_PRECONDITIONER    = "use-preconditioner";
  // const std::string ATTR_PRECONDITIONER        = "preconditioner";
  // const std::string ATTR_JACOBI_BLOCK_SIZE     = "jacobi-block-size";
  // const std::string ATTR_MAX_ITERATIONS        = "max-iterations";

  // mapping constraint
  Mapping::Constraint constraintValue{};

  mesh::PtrMeshConfiguration _meshConfig;

  // main data structure storing the configurations
  std::vector<ConfiguredMapping> _mappings;

  // Relevant information in order to store the RBF related settings,
  // as we can only instantiate the RBF classes when we know the RBF
  // which is configured in the subtag
  RBFConfiguration _rbfConfig;

  struct ExecutorConfiguration {
    enum struct Executor {
      CPU,
      CUDA,
      HIP,
      OpenMP
    };

    Executor executor = Executor::CPU;
    int      deviceId{};
    int      nThreads{};
  };

  std::unique_ptr<ExecutorConfiguration> _executorConfig;

  // Settings for the iterative solvers provided by Ginkgo
  GinkgoParameter _ginkgoParameter;

  /**
   * Configures and instantiates all mappings, which do not require
   * a subtag/ a basis function. For the RBF related mappings, this class
   * stores all relevant information, but the class is not instantiated and
   * a nullptr is returned instead. The class instantiation for the RBF
   * related mappings happens in \ref xmlTagCallback() as we need to read the
   * subtag information.
   */
  ConfiguredMapping createMapping(
      const std::string &direction,
      const std::string &type,
      const std::string &fromMeshName,
      const std::string &toMeshName,
      const std::string &geoMultiscaleType,
      const std::string &geoMultiscaleAxis,
      const double &     multiscaleRadius) const;

  /**
   * Stores additional information about the requested RBF mapping such as the
   * configured polynomial and the solver type, which is not required for all
   * the other mapping types. The information is then used later when instantiating
   * the RBF mappings in \ref xmlTagCallback().
   */
  RBFConfiguration configureRBFMapping(const std::string &type,
                                       const std::string &polynomial,
                                       bool xDead, bool yDead, bool zDead,
                                       double solverRtol,
                                       double verticesPerCluster,
                                       double relativeOverlap,
                                       bool   projectToInput) const;

  void finishRBFConfiguration();

  /// Check whether a mapping to and from the same mesh already exists
  void checkDuplicates(const ConfiguredMapping &mapping);

  /// Indicates whether the mapping here requires a basis function/ subtag,
  /// given the mapping type (e.g. nearest-neighbor).
  bool requiresBasisFunction(const std::string &mappingType) const;

  /// Given a basis function name (as a string), transforms the string into an enum of the BasisFunction
  BasisFunction parseBasisFunctions(const std::string &basisFctName) const;
};
} // namespace precice::mapping
