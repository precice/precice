#pragma once

/**
 * @file RadialBasisFctBaseMapping.hpp
 * @brief Template base class shared by all RBF-based mapping strategies.
 *
 * ## Design Rationale: Template on Basis Function
 * The `RadialBasisFctBaseMapping<RADIAL_BASIS_FUNCTION_T>` class is templated on the
 * RBF kernel type (e.g., `ThinPlateSplines`, `CompactPolynomialC2`, `Gaussian`).
 * This compile-time polymorphism means:
 * - No virtual dispatch overhead in the inner evaluation loops (critical for performance).
 * - Each instantiation is fully specialized and inlined by the compiler.
 * - The kernel's properties (compact support, positive definiteness) are available as
 *   compile-time constants, enabling `if constexpr` decisions in the solver.
 *
 * ## Dead Axes
 * In preCICE, meshes are always stored in 3D, but coupling interfaces can be 2D
 * (e.g., a flat patch in the XY-plane where z=0 for all vertices).
 * In that case, including the z-dimension in the RBF distance computation would
 * produce a degenerate (rank-deficient) system. "Dead axes" are directions where
 * all vertices share the same coordinate and can therefore be safely excluded from
 * distance computations.
 *
 * Users can also manually declare axes dead via `x-dead`, `y-dead`, `z-dead` XML attributes.
 *
 * ## Polynomial Handling
 * Many RBF interpolants require a low-order polynomial term for well-posedness.
 * The polynomial adds `1 + (number of active dimensions)` degrees of freedom:
 *   - 1 constant term
 *   - 1 linear term per active (non-dead) dimension
 *
 * ## Parallel Re-partitioning (Mesh Tagging)
 * In parallel runs, each MPI rank only holds a subset of the mesh.
 * Before building the RBF system, preCICE must ensure each rank has all the
 * remote vertices whose RBF evaluation overlaps with local vertices.
 * This is done in two rounds:
 *   1. tagMeshFirstRound()  — tags all remote vertices potentially needed.
 *   2. tagMeshSecondRound() — refines the tag set based on owned vertices.
 * See the references in each function for Benjamin Uekermann's PhD thesis.
 */

#include "impl/BasisFunctions.hpp"
#include "mapping/Mapping.hpp"
#include "mesh/Filter.hpp"
#include "precice/impl/Types.hpp"

namespace precice::mapping {

/**
 * @brief Template base class for RBF-based data mappings.
 *
 * Provides shared infrastructure for constructing RBF interpolants:
 *
 *   s(x) = Σ_i λ_i * φ(‖x − x_i‖) + p(x)
 *
 * where φ is the chosen radial basis function, λ_i are interpolation coefficients,
 * and p(x) is a low-order polynomial (optional, depending on the RBF and config).
 *
 * This base handles:
 * - Dead-axis filtering (ignoring degenerate dimensions in 2D-embedded-in-3D cases)
 * - Polynomial parameter count computation
 * - Parallel mesh partitioning (two-round vertex tagging)
 *
 * Derived classes (RadialBasisFctMapping, PetRadialBasisFctMapping, PartitionOfUnityMapping)
 * implement the actual matrix assembly, decomposition, and solve steps.
 *
 * @tparam RADIAL_BASIS_FUNCTION_T The RBF kernel type (e.g. ThinPlateSplines, CompactPolynomialC2).
 *         Must satisfy the RBF concept: hasCompactSupport(), isStrictlyPositiveDefinite(),
 *         getSupportRadius() (if compact), evaluate(double radius).
 */
template <typename RADIAL_BASIS_FUNCTION_T>
class RadialBasisFctBaseMapping : public Mapping {
public:
  /**
   * @brief Constructor.
   *
   * @param[in] constraint Specifies mapping to be consistent or conservative.
   * @param[in] dimensions Dimensionality of the meshes
   * @param[in] function Radial basis function used for mapping.
   * @param[in] xDead, yDead, zDead Deactivates mapping along an axis
   */
  RadialBasisFctBaseMapping(
      Constraint                     constraint,
      int                            dimensions,
      const RADIAL_BASIS_FUNCTION_T &function,
      std::array<bool, 3>            deadAxis,
      InitialGuessRequirement        mappingType);

  ~RadialBasisFctBaseMapping() override = default;

  // Methods, which need to be implemented in a derived class

  /// Computes the mapping coefficients from the in- and output mesh.
  void computeMapping() override = 0;

  /// Removes a computed mapping.
  void clear() override = 0;

  void tagMeshFirstRound() final override;

  void tagMeshSecondRound() final override;

protected:
  /// The instantiated RBF kernel. All distance evaluations go through this object.
  RADIAL_BASIS_FUNCTION_T _basisFunction;

  /**
   * @brief Per-dimension dead-axis flags.
   *
   * `_deadAxis[d] == true` means dimension d is "dead":
   * all vertices share the same coordinate in that direction, so it contributes
   * zero to inter-vertex distances and is excluded from the RBF computation.
   *
   * Size is always equal to getDimensions() (2 or 3).
   */
  std::vector<bool> _deadAxis;

  /**
   * @brief Computes the number of polynomial degrees of freedom.
   *
   * The polynomial part of the RBF system adds extra columns/rows to the interpolation matrix.
   * The count is: 1 (constant term) + number_of_active_dimensions (one linear term per live axis).
   *
   * For example:
   *   - 3D, no dead axes  → 1 + 3 = 4 polynomial parameters (1, x, y, z)
   *   - 2D (z dead)       → 1 + 2 = 3 polynomial parameters (1, x, y)
   *   - 1D (y, z dead)    → 1 + 1 = 2 polynomial parameters (1, x)
   *
   * @note This does NOT account for the polynomial handling mode (ON/OFF/SEPARATE);
   *       callers must gate on that separately.
   *
   * @return int the polynomial degrees of freedom
   */
  int getPolynomialParameters() const;

private:
  precice::logging::Logger _log{"mapping::RadialBasisFctBaseMapping"};

  /// converts the boolean switches to a boolean vector
  void setDeadAxis(std::array<bool, 3> deadAxis);
};

// --------------------------------------------------- HEADER IMPLEMENTATIONS

template <typename RADIAL_BASIS_FUNCTION_T>
RadialBasisFctBaseMapping<RADIAL_BASIS_FUNCTION_T>::RadialBasisFctBaseMapping(
    Constraint                     constraint,
    int                            dimensions,
    const RADIAL_BASIS_FUNCTION_T &function,
    std::array<bool, 3>            deadAxis,
    InitialGuessRequirement        mappingType)
    : Mapping(constraint, dimensions, false, mappingType),
      _basisFunction(function)
{
  // Scaled-consistent mappings require full mesh connectivity (edges/triangles)
  // to evaluate surface/volume integrals for the scaling step.
  // Other constraints (consistent, conservative) only need vertex positions.
  if (isScaledConsistent()) {
    setInputRequirement(Mapping::MeshRequirement::FULL);
    setOutputRequirement(Mapping::MeshRequirement::FULL);
  } else {
    setInputRequirement(Mapping::MeshRequirement::VERTEX);
    setOutputRequirement(Mapping::MeshRequirement::VERTEX);
  }
  setDeadAxis(deadAxis);
}

template <typename RADIAL_BASIS_FUNCTION_T>
void RadialBasisFctBaseMapping<RADIAL_BASIS_FUNCTION_T>::setDeadAxis(std::array<bool, 3> deadAxis)
{
  PRECICE_ASSERT(getDimensions() <= 3);
  PRECICE_ASSERT(_deadAxis.empty());

  // Copy only the first `dimensions` entries from the 3-element dead-axis array.
  // The third entry (z) is irrelevant and ignored for 2D meshes.
  std::copy_n(deadAxis.begin(), getDimensions(), std::back_inserter(_deadAxis));

  // Warn if the user declared z-dead on an already-2D problem (no effect).
  PRECICE_WARN_IF(getDimensions() == 2 && deadAxis[2],
                  "Setting the z-axis to dead on a 2-dimensional problem has no effect. Please remove the respective mapping's \"z-dead\" attribute.");

  // Guard: at least one axis must remain active, otherwise no distance can be computed.
  PRECICE_CHECK(std::any_of(_deadAxis.begin(), _deadAxis.end(), [](const auto &ax) { return ax == false; }), "You cannot set all axes to dead for an RBF mapping. Please remove one of the respective mapping's \"x-dead\", \"y-dead\", or \"z-dead\" attributes.");
}

template <typename RADIAL_BASIS_FUNCTION_T>
int RadialBasisFctBaseMapping<RADIAL_BASIS_FUNCTION_T>::getPolynomialParameters() const
{
  PRECICE_ASSERT(_deadAxis.size() > 0);

  // Count how many dimensions are inactive (dead).
  const int deadDimensions = std::count(_deadAxis.begin(), _deadAxis.end(), true);

  // Polynomial DOF = 1 (constant) + number_of_active_dimensions (one linear term per live axis).
  // Example: 3D with no dead axis → 1 + 3 = 4, 2D (z dead) → 1 + 2 = 3.
  return 1 + getDimensions() - deadDimensions;
}

/*
 * PARALLEL RE-PARTITIONING — TWO-ROUND MESH TAGGING
 *
 * In a parallel run, each MPI rank only holds a local portion of the coupling mesh.
 * Before building the RBF interpolation system, each rank must import ALL remote
 * vertices that its RBF stencil may need. This is done in two tagging rounds:
 *
 * Round 1 (this function): conservatively tags all remote vertices that COULD be
 *   needed, based on local bounding box + support radius.
 * Round 2 (tagMeshSecondRound): refines to only vertices within reach of OWNED vertices.
 *
 * Reference: Figure 69 in Benjamin Uekermann's PhD thesis (page 89):
 *   https://mediatum.ub.tum.de/doc/1320661/document.pdf
 */
template <typename RADIAL_BASIS_FUNCTION_T>
void RadialBasisFctBaseMapping<RADIAL_BASIS_FUNCTION_T>::tagMeshFirstRound()
{
  PRECICE_TRACE();

  // Determine which mesh is "remote" (needs filtering/tagging) and which is "local".
  // For CONSERVATIVE mappings the roles are swapped compared to CONSISTENT:
  //   consistent:  input = source (remote or local), output = target (local)
  //   conservative: input = target (local, owned by this rank), output = source (remote, needs filter)
  mesh::PtrMesh filterMesh, otherMesh;
  if (hasConstraint(CONSERVATIVE)) {
    filterMesh = output(); // remote mesh to be filtered
    otherMesh  = input();  // local mesh whose bounding box drives the filter
  } else {
    filterMesh = input();  // remote mesh to be filtered
    otherMesh  = output(); // local mesh whose bounding box drives the filter
  }

  // Ranks that have no local interface vertices don't participate in the RBF system.
  if (otherMesh->empty())
    return;

  if (_basisFunction.hasCompactSupport()) {
    // COMPACT SUPPORT: Only tag remote vertices within (local_BB + support_radius).
    // The expanded bounding box ensures that every remote vertex whose RBF kernel
    // overlaps with any local vertex is included.
    auto bb = otherMesh->getBoundingBox();
    bb.expandBy(_basisFunction.getSupportRadius());

    // Linear scan (no index tree): building an R-tree on the unfiltered global mesh
    // is more expensive than a simple loop at this stage.
    auto &vertices = filterMesh->vertices();
    std::for_each(vertices.begin(), vertices.end(), [&bb](auto &v) {
      if (bb.contains(v)) {
        v.tag(); // Mark this remote vertex as needed by the local rank
      }
    });
  } else {
    // GLOBAL SUPPORT: Every remote vertex interacts with every local vertex.
    // Must tag all remote vertices — no spatial filtering is possible.
    filterMesh->tagAll();
  }
}

/*
 * SECOND ROUND OF PARALLEL RE-PARTITIONING
 *
 * After round 1 tagged all potentially needed remote vertices, round 2 refines the set.
 * It builds a bounding box only around OWNED vertices on this rank, expands by the
 * support radius, and tags any vertex of the remote mesh that falls inside.
 *
 * This two-round process ensures the correct overlap for compact-support RBF systems
 * without over-fetching data.
 *
 * Reference: Figure 69 in Benjamin Uekermann's PhD thesis (page 89):
 *   https://mediatum.ub.tum.de/doc/1320661/document.pdf
 */
template <typename RADIAL_BASIS_FUNCTION_T>
void RadialBasisFctBaseMapping<RADIAL_BASIS_FUNCTION_T>::tagMeshSecondRound()
{
  PRECICE_TRACE();

  // Global-support RBFs: all vertices are already tagged in round 1 (tagAll), nothing to refine.
  if (not _basisFunction.hasCompactSupport())
    return;

  // Select the mesh whose vertices we want to (further) tag.
  // Same constraint-based role reversal as in tagMeshFirstRound.
  mesh::PtrMesh mesh;
  if (hasConstraint(CONSERVATIVE)) {
    mesh = output();
  } else {
    mesh = input();
  }

  // Build a bounding box that tightly covers ONLY the vertices owned by this rank.
  // Ghost/halo vertices are excluded — they may already be tagged but do not drive
  // the local RBF stencil.
  mesh::BoundingBox bb(mesh->getDimensions());
  for (mesh::Vertex &v : mesh->vertices()) {
    if (v.isOwner()) {
      PRECICE_ASSERT(v.isTagged()); // Owned vertices must have been tagged in round 1
      bb.expandBy(v);
    }
  }

  // Expand the owned bounding box by the support radius.
  // Any vertex of the remote mesh within this expanded box may contribute
  // to an RBF evaluation at one of this rank's owned vertices.
  bb.expandBy(_basisFunction.getSupportRadius());

  // Use the spatial index (R-tree) for an efficient range query rather than linear scan.
  auto vertices = mesh->index().getVerticesInsideBox(bb);
  std::for_each(vertices.begin(), vertices.end(), [&mesh](size_t v) { mesh->vertex(v).tag(); });
}
} // namespace precice::mapping
