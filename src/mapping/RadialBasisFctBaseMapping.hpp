#pragma once

#include "impl/BasisFunctions.hpp"
#include "mapping/Mapping.hpp"
#include "mesh/Filter.hpp"
#include "precice/types.hpp"
#include "utils/Event.hpp"

namespace precice {
extern bool syncMode;

namespace mapping {

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
      Constraint              constraint,
      int                     dimensions,
      RADIAL_BASIS_FUNCTION_T function,
      bool                    xDead,
      bool                    yDead,
      bool                    zDead);

  virtual ~RadialBasisFctBaseMapping() = default;

  /// Computes the mapping coefficients from the in- and output mesh.
  virtual void computeMapping() = 0;

  /// Returns true, if computeMapping() has been called.
  virtual bool hasComputedMapping() const final;

  /// Removes a computed mapping.
  virtual void clear() override;

  /// Maps input data to output data from input mesh to output mesh.
  virtual void map(int inputDataID, int outputDataID) final;

  virtual void tagMeshFirstRound() final;

  virtual void tagMeshSecondRound() final;

protected:
  /// Radial basis function type used in interpolation.
  RADIAL_BASIS_FUNCTION_T _basisFunction;

  /// true if the mapping along some axis should be ignored
  std::vector<bool> _deadAxis;

  bool _hasComputedMapping = false;
  virtual void mapConservative(int inputDataID, int outputDataID, int polyparams) = 0;
  virtual void mapConsistent(int inputDataID, int outputDataID, int polyparams)   = 0;

private:
  precice::logging::Logger _log{"mapping::RadialBasisFctBaseMapping"};


  void setDeadAxis(bool xDead, bool yDead, bool zDead);
};

// --------------------------------------------------- HEADER IMPLEMENTATIONS

template <typename RADIAL_BASIS_FUNCTION_T>
RadialBasisFctBaseMapping<RADIAL_BASIS_FUNCTION_T>::RadialBasisFctBaseMapping(
    Constraint              constraint,
    int                     dimensions,
    RADIAL_BASIS_FUNCTION_T function,
    bool                    xDead,
    bool                    yDead,
    bool                    zDead)
    : Mapping(constraint, dimensions),
      _basisFunction(function)
{
  if (constraint == SCALEDCONSISTENT) {
    setInputRequirement(Mapping::MeshRequirement::FULL);
    setOutputRequirement(Mapping::MeshRequirement::FULL);
  } else {
    setInputRequirement(Mapping::MeshRequirement::VERTEX);
    setOutputRequirement(Mapping::MeshRequirement::VERTEX);
  }
  setDeadAxis(xDead, yDead, zDead);
}

template <typename RADIAL_BASIS_FUNCTION_T>
void RadialBasisFctBaseMapping<RADIAL_BASIS_FUNCTION_T>::setDeadAxis(bool xDead, bool yDead, bool zDead)
{
  _deadAxis = {xDead, yDead};

  if (getDimensions() == 3) {
    _deadAxis.emplace_back(zDead);
  } else {
    PRECICE_ASSERT(getDimensions() == 2, "Unknown dimension.");
    if (zDead) {
      PRECICE_WARN("Setting the z-axis to dead on a 2-dimensional problem has no effect.");
    }
  }
  PRECICE_CHECK(std::any_of(_deadAxis.begin(), _deadAxis.end(), [](const auto &ax) { return ax == false; }), "You cannot choose all axes to be dead for a RBF mapping");
}

template <typename RADIAL_BASIS_FUNCTION_T>
bool RadialBasisFctBaseMapping<RADIAL_BASIS_FUNCTION_T>::hasComputedMapping() const
{
  return _hasComputedMapping;
}

template <typename RADIAL_BASIS_FUNCTION_T>
void RadialBasisFctBaseMapping<RADIAL_BASIS_FUNCTION_T>::clear()
{
  _hasComputedMapping = false;
}

template <typename RADIAL_BASIS_FUNCTION_T>
void RadialBasisFctBaseMapping<RADIAL_BASIS_FUNCTION_T>::map(
    int inputDataID,
    int outputDataID)
{
  PRECICE_TRACE(inputDataID, outputDataID);

  precice::utils::Event e("map.rbf.mapData.From" + input()->getName() + "To" + output()->getName(), precice::syncMode);

  PRECICE_ASSERT(_hasComputedMapping);
  PRECICE_ASSERT(input()->getDimensions() == output()->getDimensions(),
                 input()->getDimensions(), output()->getDimensions());
  PRECICE_ASSERT(getDimensions() == output()->getDimensions(),
                 getDimensions(), output()->getDimensions());
  {
    int valueDim = input()->data(inputDataID)->getDimensions();
    PRECICE_ASSERT(valueDim == output()->data(outputDataID)->getDimensions(),
                   valueDim, output()->data(outputDataID)->getDimensions());
  }
  int deadDimensions = 0;
  for (int d = 0; d < getDimensions(); d++) {
    if (_deadAxis[d])
      deadDimensions += 1;
  }
  int polyparams = 1 + getDimensions() - deadDimensions;

  if (hasConstraint(CONSERVATIVE)) {
    mapConservative(inputDataID, outputDataID, polyparams);
  } else {
    mapConsistent(inputDataID, outputDataID, polyparams);
  }
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

  if (otherMesh->vertices().empty())
    return; // Ranks not at the interface should never hold interface vertices

  // Tags all vertices that are inside otherMesh's bounding box, enlarged by the support radius

  if (_basisFunction.hasCompactSupport()) {
    auto bb = otherMesh->getBoundingBox();
    bb.expandBy(_basisFunction.getSupportRadius());

    auto vertices = filterMesh->index().getVerticesInsideBox(bb);
    std::for_each(vertices.begin(), vertices.end(), [&filterMesh](size_t v) { filterMesh->vertices()[v].tag(); });
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
  std::for_each(vertices.begin(), vertices.end(), [&mesh](size_t v) { mesh->vertices()[v].tag(); });
}
} // namespace mapping
} // namespace precice
