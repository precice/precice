#pragma once

#include <string>
#include <vector>
#include "logging/Logger.hpp"
#include "mapping/SharedPointer.hpp"
#include "mesh/SharedPointer.hpp"
#include "xml/XMLTag.hpp"

namespace precice {
namespace mapping {

/// How to handle the polynomial?
/**
 * ON: Include it in the system matrix
 * OFF: Omit it altogether
 * SEPARATE: Compute it separately using least-squares QR.
 */
enum class Polynomial {
  ON,
  OFF,
  SEPARATE
};

enum class Preallocation {
  OFF,
  COMPUTE,
  ESTIMATE,
  SAVE,
  TREE
};

enum class RBFType {
  EIGEN,
  PETSc
};

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
    bool isRBF;
  };

  struct RBFParameter {

    enum struct Type {
      ShapeParameter,
      SupportRadius
    };

    Type   type{};
    double value{};
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

  void resetMappings()
  {
    _mappings.clear();
  }

private:
  mutable logging::Logger _log{"config:MappingConfiguration"};

  const std::string TAG = "mapping";

  const std::string ATTR_DIRECTION      = "direction";
  const std::string ATTR_FROM           = "from";
  const std::string ATTR_TO             = "to";
  const std::string ATTR_TYPE           = "type";
  const std::string ATTR_CONSTRAINT     = "constraint";
  const std::string ATTR_SHAPE_PARAM    = "shape-parameter";
  const std::string ATTR_SUPPORT_RADIUS = "support-radius";
  const std::string ATTR_SOLVER_RTOL    = "solver-rtol";
  const std::string ATTR_X_DEAD         = "x-dead";
  const std::string ATTR_Y_DEAD         = "y-dead";
  const std::string ATTR_Z_DEAD         = "z-dead";
  const std::string ATTR_USE_QR         = "use-qr-decomposition";

  const std::string ATTR_VERTICES_PER_PARTITION = "vertices-per-partition";
  const std::string ATTR_RELATIVE_OVERLAP       = "relative-overlap";
  const std::string ATTR_PROJECT_TO_INPUT       = "project-to-input";

  const std::string VALUE_WRITE                     = "write";
  const std::string VALUE_READ                      = "read";
  const std::string VALUE_CONSISTENT                = "consistent";
  const std::string VALUE_CONSERVATIVE              = "conservative";
  const std::string VALUE_SCALED_CONSISTENT_SURFACE = "scaled-consistent-surface";
  const std::string VALUE_SCALED_CONSISTENT_VOLUME  = "scaled-consistent-volume";

  const std::string VALUE_NEAREST_NEIGHBOR          = "nearest-neighbor";
  const std::string VALUE_NEAREST_PROJECTION        = "nearest-projection";
  const std::string VALUE_LINEAR_CELL_INTERPOLATION = "linear-cell-interpolation";

  const std::string VALUE_RBF_TPS               = "rbf-thin-plate-splines";
  const std::string VALUE_RBF_MULTIQUADRICS     = "rbf-multiquadrics";
  const std::string VALUE_RBF_INV_MULTIQUADRICS = "rbf-inverse-multiquadrics";
  const std::string VALUE_RBF_VOLUME_SPLINES    = "rbf-volume-splines";
  const std::string VALUE_RBF_GAUSSIAN          = "rbf-gaussian";
  const std::string VALUE_RBF_CTPS_C2           = "rbf-compact-tps-c2";
  const std::string VALUE_RBF_CPOLYNOMIAL_C0    = "rbf-compact-polynomial-c0";
  const std::string VALUE_RBF_CPOLYNOMIAL_C2    = "rbf-compact-polynomial-c2";
  const std::string VALUE_RBF_CPOLYNOMIAL_C4    = "rbf-compact-polynomial-c4";
  const std::string VALUE_RBF_CPOLYNOMIAL_C6    = "rbf-compact-polynomial-c6";
  const std::string VALUE_RBF_POU               = "rbf-pou";

  const std::string VALUE_NEAREST_NEIGHBOR_GRADIENT = "nearest-neighbor-gradient";

  mesh::PtrMeshConfiguration _meshConfig;

  std::vector<ConfiguredMapping> _mappings;

  ConfiguredMapping createMapping(
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
      Preallocation                    preallocation,
      double                           verticesPerPartition,
      double                           relativeOverlap,
      bool                             projectToInput) const;

  /// Check whether a mapping to and from the same mesh already exists
  void checkDuplicates(const ConfiguredMapping &mapping);
};
} // namespace mapping
} // namespace precice
