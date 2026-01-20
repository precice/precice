#pragma once

#include "impl/BasisFunctions.hpp"
#include "mapping/Mapping.hpp"
#include "mesh/Filter.hpp"
#include "precice/impl/Types.hpp"

namespace precice::mapping {

/**
 * @brief Mapping with radial basis functions.
 *
 * With help of the input data points and values an interpolant is constructed.
 * The interpolant is formed by a weighted sum of conditionally positive radial
 * basis functions and a (low order) polynomial, and evaluated at the output
 * data points.
 *
 * The radial basis function type has to be given as template parameter, and has
 * to be one of the defined types in this file.
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
  /// Radial basis function type used in interpolation.
  RADIAL_BASIS_FUNCTION_T _basisFunction;

  /// true if the mapping along some axis should be ignored
  std::vector<bool> _deadAxis;

  /**
   * @brief Computes the number of polynomial degrees of freedom based on the problem dimension and the dead axis
   *
   * @note This function does not take the handling of the polynomial (ON, OFF, SEPARATE) into account.
   *
   * @return int the polynomial degrees of freedom of the RBF system
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
  std::copy_n(deadAxis.begin(), getDimensions(), std::back_inserter(_deadAxis));
  PRECICE_WARN_IF(getDimensions() == 2 && deadAxis[2],
                  "Setting the z-axis to dead on a 2-dimensional problem has no effect. Please remove the respective mapping's \"z-dead\" attribute.");
  PRECICE_CHECK(std::any_of(_deadAxis.begin(), _deadAxis.end(), [](const auto &ax) { return ax == false; }), "You cannot set all axes to dead for an RBF mapping. Please remove one of the respective mapping's \"x-dead\", \"y-dead\", or \"z-dead\" attributes.");
}

template <typename RADIAL_BASIS_FUNCTION_T>
int RadialBasisFctBaseMapping<RADIAL_BASIS_FUNCTION_T>::getPolynomialParameters() const
{
  PRECICE_ASSERT(_deadAxis.size() > 0);
  // Count the dead axis
  const int deadDimensions = std::count(_deadAxis.begin(), _deadAxis.end(), true);
  // Formula for the polynomial parameters
  return 1 + getDimensions() - deadDimensions;
}

/*
 * For the re-partitioning process with RBF mappings, also compare Figure 69 in Benjamin U's thesis (page 89).
 * https://mediatum.ub.tum.de/doc/1320661/document.pdf
 */
template <typename RADIAL_BASIS_FUNCTION_T>
void RadialBasisFctBaseMapping<RADIAL_BASIS_FUNCTION_T>::tagMeshFirstRound()
{
  PRECICE_TRACE();
  mesh::PtrMesh filterMesh, otherMesh;
  if (hasConstraint(CONSERVATIVE)) {
    filterMesh = output(); // remote
    otherMesh  = input();  // local
  } else {
    filterMesh = input();  // remote
    otherMesh  = output(); // local
  }

  if (otherMesh->empty())
    return; // Ranks not at the interface should never hold interface vertices

  // Tags all vertices that are inside otherMesh's bounding box, enlarged by the support radius
  if (_basisFunction.hasCompactSupport()) {
    auto bb = otherMesh->getBoundingBox();
    bb.expandBy(_basisFunction.getSupportRadius());

    // We don't make use of the index tree here, because constructing the index tree on the
    // (unfiltered) mesh is expensive
    auto &vertices = filterMesh->vertices();
    std::for_each(vertices.begin(), vertices.end(), [&bb](auto &v) {
      if (bb.contains(v)) {
        v.tag();
      }
    });
  } else {
    filterMesh->tagAll();
  }
}

/*
 * For the re-partitioning process with RBF mappings, also compare Figure 69 in Benjamin U's thesis (page 89).
 * https://mediatum.ub.tum.de/doc/1320661/document.pdf
 */
template <typename RADIAL_BASIS_FUNCTION_T>
void RadialBasisFctBaseMapping<RADIAL_BASIS_FUNCTION_T>::tagMeshSecondRound()
{
  PRECICE_TRACE();

  if (not _basisFunction.hasCompactSupport())
    return; // Tags should not be changed

  mesh::PtrMesh mesh; // The mesh we want to filter

  if (hasConstraint(CONSERVATIVE)) {
    mesh = output();
  } else {
    mesh = input();
  }

  mesh::BoundingBox bb(mesh->getDimensions());

  // Construct bounding box around all owned vertices
  for (mesh::Vertex &v : mesh->vertices()) {
    if (v.isOwner()) {
      PRECICE_ASSERT(v.isTagged()); // Should be tagged from the first round
      bb.expandBy(v);
    }
  }
  // Enlarge bb by support radius
  bb.expandBy(_basisFunction.getSupportRadius());
  auto vertices = mesh->index().getVerticesInsideBox(bb);
  std::for_each(vertices.begin(), vertices.end(), [&mesh](size_t v) { mesh->vertex(v).tag(); });
}
} // namespace precice::mapping
