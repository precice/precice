#include "RadialGeoMultiscaleMapping.hpp"
#include <Eigen/src/Core/Matrix.h>
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
  PRECICE_TRACE(output()->vertices().size());

  PRECICE_ASSERT(input().get() != nullptr);
  PRECICE_ASSERT(output().get() != nullptr);

  size_t const inSize  = input()->vertices().size();
  size_t const outSize = output()->vertices().size();

  int effectiveCoordinate = static_cast<std::underlying_type_t<MultiscaleType>>(_axis);
  PRECICE_ASSERT(effectiveCoordinate == static_cast<std::underlying_type_t<MultiscaleType>>(MultiscaleAxis::X) ||
                     effectiveCoordinate == static_cast<std::underlying_type_t<MultiscaleType>>(MultiscaleAxis::Y) ||
                     effectiveCoordinate == static_cast<std::underlying_type_t<MultiscaleType>>(MultiscaleAxis::Z),
                 "Unknown multiscale axis type.")

  if (getConstraint() == CONSISTENT) {
    PRECICE_DEBUG("Compute consistent mapping");
    if (_type == MultiscaleType::SPREAD) {
      // Determine principle axis midpoints as borders to assign 3D vertices to their respective 1D vertex
      Eigen::VectorXd axisMidpoints(inSize);
      auto &          inputVerticesRef = input()->vertices();
      // Order the vertices of the 1D input mesh, to correctly compute the midpoints of the segments between neighboring vertices
      std::vector<size_t> ordered_vertex_indices(input()->vertices().size());
      std::iota(ordered_vertex_indices.begin(), ordered_vertex_indices.end(), 0);
      std::stable_sort(ordered_vertex_indices.begin(), ordered_vertex_indices.end(),
                       [effectiveCoordinate, &inputVerticesRef](const size_t aindex, const size_t bindex) {
                         return inputVerticesRef[aindex].rawCoords()[effectiveCoordinate] < inputVerticesRef[bindex].rawCoords()[effectiveCoordinate];
                       });
      // Compute the midpoints of the 1D input mesh
      for (size_t i = 0; i < (inSize - 1); i++) {
        auto axisPositionCurrent = inputVerticesRef[ordered_vertex_indices[i]].rawCoords()[effectiveCoordinate];
        auto axisPositionNext    = inputVerticesRef[ordered_vertex_indices[i] + 1].rawCoords()[effectiveCoordinate];
        axisMidpoints(i)         = (axisPositionCurrent + axisPositionNext) / 2;
      }
      axisMidpoints(inSize - 1) = std::numeric_limits<double>::max(); // large number, such that vertices after the last midpoint are still assigned

      /*
        3D vertices are projected onto the 1D axis and the data is then mapped
        to the nearest neighbors of the 1D vertices in projection space.
      */
      _vertexIndicesSpread.clear();
      _vertexIndicesSpread.reserve(output()->vertices().size());
      for (size_t i = 0; i < outSize; i++) {
        auto vertexCoord = output()->vertices()[i].rawCoords()[effectiveCoordinate];
        int  index       = 0;
        while (vertexCoord > axisMidpoints(index)) {
          PRECICE_ASSERT(index + 1 < static_cast<int>(inSize));
          ++index;
        }
        _vertexIndicesSpread.push_back(index);
      }
    } else {
      PRECICE_ASSERT(_type == MultiscaleType::COLLECT)
      /*
        3D vertices are projected onto the 1D axis and the data is then mapped
        to (and averaged at) the nearest 1D vertex in projection space.
      */

      // determine principle axis midpoints as borders to assign 3D vertices to their respective 1D vertex
      Eigen::VectorXd axisMidpoints(outSize);
      auto &          outputVerticesRef = output()->vertices();

      // Order the vertices of the 1D output mesh, to correctly compute the midpoints of the segments between neighboring vertices
      std::vector<size_t> ordered_vertex_indices(output()->vertices().size());
      std::iota(ordered_vertex_indices.begin(), ordered_vertex_indices.end(), 0);
      std::stable_sort(ordered_vertex_indices.begin(), ordered_vertex_indices.end(),
                       [effectiveCoordinate, &outputVerticesRef](const size_t aindex, const size_t bindex) {
                         return outputVerticesRef[aindex].rawCoords()[effectiveCoordinate] < outputVerticesRef[bindex].rawCoords()[effectiveCoordinate];
                       });
      // Compute the midpoints of the 1D output mesh
      for (size_t i = 0; i < (outSize - 1); i++) {
        auto axisPositionCurrent = output()->vertices()[i].rawCoords()[effectiveCoordinate];
        auto axisPositionNext    = output()->vertices()[i + 1].rawCoords()[effectiveCoordinate];
        axisMidpoints(i)         = (axisPositionCurrent + axisPositionNext) / 2;
      }
      axisMidpoints(outSize - 1) = std::numeric_limits<double>::max(); // large number, such that vertices after the last midpoint are still assigned

      std::vector<int> counters(outSize); // counts number of vertices in between midpoints for averaging

      // Identify which vertex (index) of the 3D mesh corresponds to which vertex (index) of the 1D meesh
      _vertexIndicesCollect.clear();
      _vertexIndicesCollect.reserve(input()->vertices().size());
      for (size_t i = 0; i < inSize; i++) {
        auto vertexCoords = input()->vertices()[i].rawCoords()[effectiveCoordinate];
        int  index        = 0;
        while (vertexCoords > axisMidpoints(index)) {
          PRECICE_ASSERT(index + 1 < static_cast<int>(outSize));
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
  Eigen::VectorXd &      outputValues     = outData;

  size_t const inSize  = input()->vertices().size();
  size_t const outSize = output()->vertices().size();

  PRECICE_ASSERT(!output()->vertices().empty());
  auto outDataDimensions = outputValues.size() / outSize;

  // Check that the number of values for the input and output is right according to their dimensions
  PRECICE_ASSERT((inputValues.size() / inDataDimensions == static_cast<int>(input()->vertices().size())),
                 inputValues.size(), inDataDimensions, input()->vertices().size());
  PRECICE_ASSERT((outputValues.size() / outDataDimensions == static_cast<int>(output()->vertices().size())),
                 outputValues.size(), outDataDimensions, output()->vertices().size());

  // We currently don't support 1D data, so we need that the user specifies data of the same dimensions on both sides
  PRECICE_ASSERT(inDataDimensions == outDataDimensions);

  PRECICE_DEBUG("Map consistent");
  if (_type == MultiscaleType::SPREAD) {
    // assign 1D vertex value to all 3D vertices in vicinity
    for (size_t i = 0; i < outSize; i++) {
      outputValues((i * outDataDimensions)) = inputValues(_vertexIndicesSpread[i] * inDataDimensions);
    }
  } else {
    PRECICE_ASSERT(_type == MultiscaleType::COLLECT);
    /*
      3D vertices are projected onto the 1D axis and the data is then mapped
      to (and averaged at) the nearest 1D vertex in projection space.
    */

    for (size_t i = 0; i < outSize; i++) {
      outputValues((i * outDataDimensions)) = 0;
    }
    // assign the 1D vertex the average of all 3D vertex values in vicinity
    for (size_t i = 0; i < inSize; i++) {
      outputValues(_vertexIndicesCollect[i] * outDataDimensions) += inputValues(i * inDataDimensions);
    }
    for (size_t i = 0; i < outSize; i++) {
      outputValues(i) = outputValues(i) / _vertexCounter[i];
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
