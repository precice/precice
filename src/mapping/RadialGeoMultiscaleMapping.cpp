#include "RadialGeoMultiscaleMapping.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <numeric>
#include "mesh/Mesh.hpp"

namespace precice::mapping {

RadialGeoMultiscaleMapping::RadialGeoMultiscaleMapping(
    Constraint     constraint,
    int            dimensions,
    MultiscaleType type,
    MultiscaleAxis axis)
    : Mapping(constraint, dimensions, false, Mapping::InitialGuessRequirement::None),
      _type(type), _axis(axis)
{
  setInputRequirement(Mapping::MeshRequirement::VERTEX);
  setOutputRequirement(Mapping::MeshRequirement::VERTEX);
}

void RadialGeoMultiscaleMapping::computeMapping()
{
  PRECICE_TRACE(output()->nVertices());

  PRECICE_ASSERT(input().get() != nullptr);
  PRECICE_ASSERT(output().get() != nullptr);

  size_t const inSize  = input()->nVertices();
  size_t const outSize = output()->nVertices();

  int effectiveCoordinate = static_cast<std::underlying_type_t<MultiscaleType>>(_axis);
  PRECICE_ASSERT(effectiveCoordinate == static_cast<std::underlying_type_t<MultiscaleType>>(MultiscaleAxis::X) ||
                     effectiveCoordinate == static_cast<std::underlying_type_t<MultiscaleType>>(MultiscaleAxis::Y) ||
                     effectiveCoordinate == static_cast<std::underlying_type_t<MultiscaleType>>(MultiscaleAxis::Z),
                 "Unknown multiscale axis type.");

  if (getConstraint() == CONSISTENT) {
    PRECICE_DEBUG("Compute consistent mapping");
    if (_type == MultiscaleType::SPREAD) {
      // Determine principle axis midpoints as borders to assign 3D vertices to their respective 1D vertex
      Eigen::VectorXd axisMidpoints(inSize);
      auto           &inputVerticesRef = input()->vertices();
      // Order the vertices of the 1D input mesh, to correctly compute the midpoints of the segments between neighboring vertices
      std::vector<size_t> ordered_vertex_indices(input()->nVertices());
      std::iota(ordered_vertex_indices.begin(), ordered_vertex_indices.end(), 0);
      std::stable_sort(ordered_vertex_indices.begin(), ordered_vertex_indices.end(),
                       [effectiveCoordinate, &inputVerticesRef](const size_t aindex, const size_t bindex) {
                         return inputVerticesRef[aindex].coord(effectiveCoordinate) < inputVerticesRef[bindex].coord(effectiveCoordinate);
                       });
      // Compute the midpoints of the 1D input mesh
      for (size_t i = 0; i < (inSize - 1); i++) {
        auto axisPositionCurrent = inputVerticesRef[ordered_vertex_indices[i]].coord(effectiveCoordinate);
        auto axisPositionNext    = inputVerticesRef[ordered_vertex_indices[i] + 1].coord(effectiveCoordinate);
        axisMidpoints(i)         = (axisPositionCurrent + axisPositionNext) / 2;
      }
      axisMidpoints(inSize - 1) = std::numeric_limits<double>::max(); // large number, such that vertices after the last midpoint are still assigned

      /*
        3D vertices are projected onto the 1D axis and the data is then mapped
        to the nearest neighbors of the 1D vertices in projection space.
      */
      _vertexIndicesSpread.clear();
      _vertexIndicesSpread.reserve(output()->nVertices());
      for (size_t i = 0; i < outSize; i++) {
        auto   vertexCoord = output()->vertex(i).coord(effectiveCoordinate);
        size_t index       = 0;
        while (vertexCoord > axisMidpoints(index)) {
          PRECICE_ASSERT(index + 1 < inSize);
          ++index;
        }
        _vertexIndicesSpread.push_back(index);
      }
    } else {
      PRECICE_ASSERT(_type == MultiscaleType::COLLECT);
      /*
        3D vertices are projected onto the 1D axis and the data is then mapped
        to (and averaged at) the nearest 1D vertex in projection space.
      */

      // determine principle axis midpoints as borders to assign 3D vertices to their respective 1D vertex
      Eigen::VectorXd axisMidpoints(outSize);
      auto           &outputVerticesRef = output()->vertices();

      // Order the vertices of the 1D output mesh, to correctly compute the midpoints of the segments between neighboring vertices
      std::vector<size_t> ordered_vertex_indices(output()->nVertices());
      std::iota(ordered_vertex_indices.begin(), ordered_vertex_indices.end(), 0);
      std::stable_sort(ordered_vertex_indices.begin(), ordered_vertex_indices.end(),
                       [effectiveCoordinate, &outputVerticesRef](const size_t aindex, const size_t bindex) {
                         return outputVerticesRef[aindex].coord(effectiveCoordinate) < outputVerticesRef[bindex].coord(effectiveCoordinate);
                       });
      // Compute the midpoints of the 1D output mesh
      for (size_t i = 0; i < (outSize - 1); i++) {
        auto axisPositionCurrent = output()->vertex(i).coord(effectiveCoordinate);
        auto axisPositionNext    = output()->vertex(i + 1).coord(effectiveCoordinate);
        axisMidpoints(i)         = (axisPositionCurrent + axisPositionNext) / 2;
      }
      axisMidpoints(outSize - 1) = std::numeric_limits<double>::max(); // large number, such that vertices after the last midpoint are still assigned

      std::vector<size_t> counters(outSize); // counts number of vertices in between midpoints for averaging

      // Identify which vertex (index) of the 3D mesh corresponds to which vertex (index) of the 1D meesh
      _vertexIndicesCollect.clear();
      _vertexIndicesCollect.reserve(input()->nVertices());
      for (size_t i = 0; i < inSize; i++) {
        auto   vertexCoords = input()->vertex(i).coord(effectiveCoordinate);
        size_t index        = 0;
        while (vertexCoords > axisMidpoints(index)) {
          PRECICE_ASSERT(index + 1 < outSize);
          ++index;
        }
        _vertexIndicesCollect.push_back(index);
        counters[index] += 1;
      }
      _vertexCounter = std::move(counters);
    }

  } else {
    PRECICE_ASSERT(getConstraint() == CONSERVATIVE);
    PRECICE_ASSERT(false, "Not yet implemented");
    PRECICE_DEBUG("Compute conservative mapping");
  }
  _hasComputedMapping = true;
}

void RadialGeoMultiscaleMapping::clear()
{
  PRECICE_TRACE();
  _vertexIndicesSpread.clear();
  _vertexIndicesCollect.clear();
  _hasComputedMapping = false;
}

void RadialGeoMultiscaleMapping::mapConservative(const time::Sample &inData, Eigen::VectorXd &outData)
{
  PRECICE_ASSERT(getConstraint() == CONSERVATIVE);
  PRECICE_ASSERT(false, "Not yet implemented");
  PRECICE_DEBUG("Map conservative");
}

void RadialGeoMultiscaleMapping::mapConsistent(const time::Sample &inData, Eigen::VectorXd &outData)
{
  PRECICE_TRACE();

  const int              inDataDimensions = inData.dataDims;
  const Eigen::VectorXd &inputValues      = inData.values;
  Eigen::VectorXd       &outputValues     = outData;

  size_t const inSize  = input()->nVertices();
  size_t const outSize = output()->nVertices();

  PRECICE_ASSERT(!output()->empty());
  auto outDataDimensions = outputValues.size() / outSize;

  // Check that the number of values for the input and output is right according to their dimensions
  PRECICE_ASSERT((inputValues.size() / static_cast<std::size_t>(inDataDimensions) == input()->nVertices()),
                 inputValues.size(), inDataDimensions, input()->nVertices());
  PRECICE_ASSERT((outputValues.size() / outDataDimensions == output()->nVertices()),
                 outputValues.size(), outDataDimensions, output()->nVertices());

  // We currently don't support 1D data, so we need that the user specifies data of the same dimensions on both sides
  PRECICE_ASSERT(static_cast<std::size_t>(inDataDimensions) == outDataDimensions);

  PRECICE_DEBUG("Map consistent");
  if (_type == MultiscaleType::SPREAD) {
    // assign 1D vertex value to all 3D vertices in vicinity
    for (size_t i = 0; i < outSize; i++) {
      outputValues(i * outDataDimensions) = inputValues(_vertexIndicesSpread[i] * inDataDimensions);
    }
  } else {
    PRECICE_ASSERT(_type == MultiscaleType::COLLECT);
    /*
      3D vertices are projected onto the 1D axis and the data is then mapped
      to (and averaged at) the nearest 1D vertex in projection space.
    */

    for (size_t i = 0; i < outSize; i++) {
      outputValues(i * outDataDimensions) = 0;
    }
    // assign the 1D vertex the average of all 3D vertex values in vicinity
    for (size_t i = 0; i < inSize; i++) {
      outputValues(_vertexIndicesCollect[i] * outDataDimensions) += inputValues(i * inDataDimensions);
    }
    for (size_t i = 0; i < outSize; i++) {
      outputValues(i * outDataDimensions) = outputValues(i * outDataDimensions) / _vertexCounter[i];
    }
  }
}

void RadialGeoMultiscaleMapping::tagMeshFirstRound()
{
  PRECICE_TRACE();

  computeMapping();

  if (getConstraint() == CONSISTENT) {
    PRECICE_ASSERT(_type == MultiscaleType::SPREAD, "Not yet implemented");

    // tag all vertices of the 1D mesh
    input()->tagAll();

  } else {
    PRECICE_ASSERT(getConstraint() == CONSERVATIVE);
    PRECICE_ASSERT(false, "Not yet implemented");
  }

  clear();
}

void RadialGeoMultiscaleMapping::tagMeshSecondRound()
{
  PRECICE_TRACE();
  // no operation needed here for the moment
}

std::string RadialGeoMultiscaleMapping::getName() const
{
  return "radial-geomultiscale";
}

} // namespace precice::mapping
