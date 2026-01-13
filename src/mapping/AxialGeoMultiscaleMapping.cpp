#include "AxialGeoMultiscaleMapping.hpp"
#include <Eigen/Core>
#include "mesh/Mesh.hpp"

namespace precice::mapping {

AxialGeoMultiscaleMapping::AxialGeoMultiscaleMapping(
    Constraint             constraint,
    int                    dimensions,
    MultiscaleDimension    dimension,
    MultiscaleType         type,
    MultiscaleAxis         axis,
    double                 radius,
    SpreadProfile          profile,
    MultiscaleCrossSection crossSection)
    : Mapping(constraint, dimensions, false, Mapping::InitialGuessRequirement::None),
      _dimension(dimension),
      _type(type),
      _axis(axis),
      _radius(radius),
      _profile(profile),
      _crossSection(crossSection)
{
  setInputRequirement(Mapping::MeshRequirement::VERTEX);
  setOutputRequirement(Mapping::MeshRequirement::VERTEX);
}

void AxialGeoMultiscaleMapping::computeMapping()
{
  PRECICE_TRACE(output()->nVertices());

  PRECICE_ASSERT(input().get() != nullptr);
  PRECICE_ASSERT(output().get() != nullptr);

  if (getConstraint() == CONSISTENT) {
    PRECICE_DEBUG("Compute consistent mapping");
    if (_type == MultiscaleType::SPREAD) {
      const int outDataDimensions   = 3;
      int       effectiveCoordinate = static_cast<std::underlying_type_t<MultiscaleType>>(_axis); // Convert enum struct to int
      PRECICE_ASSERT(effectiveCoordinate == static_cast<std::underlying_type_t<MultiscaleType>>(MultiscaleAxis::X) ||
                         effectiveCoordinate == static_cast<std::underlying_type_t<MultiscaleType>>(MultiscaleAxis::Y) ||
                         effectiveCoordinate == static_cast<std::underlying_type_t<MultiscaleType>>(MultiscaleAxis::Z),
                     "Unknown multiscale axis type.");
      size_t const inSize  = input()->nVertices();
      size_t const outSize = output()->nVertices();
      if (_dimension == MultiscaleDimension::D1D3 || _dimension == MultiscaleDimension::D1D2) {
        PRECICE_CHECK(input()->nVertices() == 1, "You can only define an axial geometric multiscale 1D-{} mapping of type spread from a mesh with exactly one vertex.");

        /* When we add support for 1D meshes (https://github.com/precice/precice/issues/1669),
          we should check for valid dimensions combination.

          const int inDataDimensions  = input()->getDimensions();
          const int outDataDimensions = output()->getDimensions();

          PRECICE_CHECK(input()->getDimensions() == 1, "The input mesh on an axial geometric multiscale mapping can only be 1D at the moment, but it was defined to be {}.", input()->getDimensions())
          PRECICE_CHECK(output()->getDimensions() == 3, "The output mesh on an axial geometric multiscale mapping can only be 3D at the moment, but it was defined to be {}.", input()->getDimensions())
        */

        // compute distances between 1D vertex and 3D vertices
        mesh::Vertex    &v0                           = input()->vertex(0);
        constexpr double distance_to_radius_threshold = 1.05;
        if (_crossSection == MultiscaleCrossSection::CIRCLE) {
          _vertexDistances.clear();
          _vertexDistances.reserve(output()->nVertices());

          for (size_t i = 0; i < outSize; i++) {
            Eigen::VectorXd difference(outDataDimensions);
            difference = v0.getCoords();
            difference -= output()->vertex(i).getCoords();
            double distance_to_radius = difference.norm() / _radius;
            PRECICE_CHECK(distance_to_radius <= distance_to_radius_threshold, "Output mesh has vertices that do not coincide with the geometric multiscale interface defined by the input mesh. Ratio of vertex distance to radius is {} (which is larger than the assumed threshold of distance_to_radius_threshold).", distance_to_radius);
            _vertexDistances.push_back(distance_to_radius);
          }
        } else {
          PRECICE_ASSERT(_crossSection == MultiscaleCrossSection::SQUARE);
          _vertexTransverseCoords.clear();
          _vertexTransverseCoords.reserve(output()->nVertices());

          const int t1 = (effectiveCoordinate + 1) % 3;
          const int t2 = (effectiveCoordinate + 2) % 3;

          for (size_t i = 0; i < outSize; i++) {
            Eigen::VectorXd difference(outDataDimensions);
            difference = v0.getCoords();
            difference -= output()->vertex(i).getCoords();

            const double s1 = difference[t1] / _radius;
            const double s2 = difference[t2] / _radius;

            const double squareNorm = std::max(std::abs(s1), std::abs(s2));
            PRECICE_CHECK(squareNorm <= distance_to_radius_threshold,
                          "Output mesh has vertices that do not coincide with the square multiscale interface defined by the input mesh. "
                          "max(|xn|,|yn|) is {} (threshold {}).",
                          squareNorm, distance_to_radius_threshold);
            _vertexTransverseCoords.push_back({s1, s2});
          }
        }

      } else {
        PRECICE_ASSERT(_dimension == MultiscaleDimension::D2D3);
        PRECICE_CHECK(input()->nVertices() > 1, "You can only define an axial geometric multiscale 2D-3D mapping of type spread from a mesh with more than 1 vertex.");
        Eigen::Vector3d minC = input()->vertex(0).getCoords().head<3>();
        Eigen::Vector3d maxC = minC;
        for (size_t i = 1; i < input()->nVertices(); ++i) {
          const Eigen::Vector3d c = input()->vertex(i).getCoords().head<3>();
          minC                    = minC.cwiseMin(c);
          maxC                    = maxC.cwiseMax(c);
        }
        const Eigen::Vector3d span = maxC - minC;
        _lineCoord                 = 0;
        if (span[1] > span[_lineCoord])
          _lineCoord = 1;
        if (span[2] > span[_lineCoord])
          _lineCoord = 2;
        _nearestVertex.clear();
        _nearestVertex.reserve(output()->nVertices());
        _vertexDistances.clear();
        _vertexDistances.reserve(output()->nVertices());
        _maxDistancePerInput.clear();
        _maxDistancePerInput.resize(input()->nVertices(), 0.0);
        for (size_t j = 0; j < outSize; j++) {
          const Eigen::VectorXd &xOut          = output()->vertex(j).getCoords();
          double                 bestDistance2 = std::numeric_limits<double>::max();
          int                    bestIdx       = -1;
          for (size_t i = 0; i < inSize; i++) {
            const Eigen::VectorXd &xIn        = input()->vertex(i).getCoords();
            Eigen::VectorXd        difference = xOut - xIn;
            double                 distance2  = difference.squaredNorm();
            if (distance2 < bestDistance2) {
              bestDistance2 = distance2;
              bestIdx       = static_cast<int>(i);
            }
          }
          PRECICE_ASSERT(bestIdx >= 0, "Could not find nearest input vertex for output vertex {}", j);
          _nearestVertex.push_back(bestIdx);
          double distance = std::sqrt(bestDistance2);
          _vertexDistances.push_back(distance);
          if (distance > _maxDistancePerInput[static_cast<size_t>(bestIdx)]) {
            _maxDistancePerInput[static_cast<size_t>(bestIdx)] = distance;
          }
        }
      }
    } else {
      PRECICE_ASSERT(_type == MultiscaleType::COLLECT);
      size_t const inSize  = input()->nVertices();
      size_t const outSize = output()->nVertices();
      if (_dimension == MultiscaleDimension::D1D3) {
        PRECICE_CHECK(output()->nVertices() == 1, "You can only define an axial geometric multiscale 1D-{} mapping of type collect to a mesh with exactly one vertex.");
      } else if (_dimension == MultiscaleDimension::D1D2) {
        PRECICE_CHECK(output()->nVertices() == 1,
                      "You can only define an axial geometric multiscale 1D-2D mapping of type collect to a mesh with exactly one vertex.");
      } else {
        PRECICE_ASSERT(_dimension == MultiscaleDimension::D2D3);
        PRECICE_CHECK(outSize > 1, "You can only define an axial geometric multiscale 2D-3D mapping of type collect to a mesh with more than 1 vertex.");
        Eigen::Vector3d minC = output()->vertex(0).getCoords().head<3>();
        Eigen::Vector3d maxC = minC;
        for (size_t i = 1; i < output()->nVertices(); ++i) {
          const Eigen::Vector3d c = output()->vertex(i).getCoords().head<3>();
          minC                    = minC.cwiseMin(c);
          maxC                    = maxC.cwiseMax(c);
        }
        const Eigen::Vector3d span = maxC - minC;
        _lineCoord                 = 0;
        if (span[1] > span[_lineCoord])
          _lineCoord = 1;
        if (span[2] > span[_lineCoord])
          _lineCoord = 2;
        _collectBands.clear();
        _collectBands.resize(output()->nVertices());
        for (size_t i = 0; i < inSize; i++) {
          const Eigen::VectorXd &xIn          = input()->vertex(i).getCoords();
          double                 bestDistance = std::numeric_limits<double>::max();
          int                    bestIdx      = -1;
          for (size_t j = 0; j < outSize; j++) {
            const Eigen::VectorXd &xOut       = output()->vertex(j).getCoords();
            Eigen::VectorXd        difference = xIn - xOut;
            double                 distance   = difference.squaredNorm();
            if (distance < bestDistance) {
              bestDistance = distance;
              bestIdx      = static_cast<int>(j);
            }
          }
          PRECICE_ASSERT(bestIdx >= 0, "Could not find nearest output vertex for input vertex {}", i);
          _collectBands[bestIdx].push_back(static_cast<int>(i));
        }
      }
    }

  } else {
    PRECICE_ASSERT(getConstraint() == CONSERVATIVE);
    PRECICE_UNREACHABLE("Axial conservative geometric multiscale mapping is not implemented");
    PRECICE_DEBUG("Compute conservative mapping");
  }
  _hasComputedMapping = true;
}

void AxialGeoMultiscaleMapping::clear()
{
  PRECICE_TRACE();
  _vertexDistances.clear();
  _hasComputedMapping = false;
}

void AxialGeoMultiscaleMapping::mapConservative(const time::Sample &inData, Eigen::VectorXd &outData)
{
  PRECICE_ASSERT(getConstraint() == CONSERVATIVE);
  PRECICE_UNREACHABLE("Axial conservative geometric multiscale mapping is not implemented");
  PRECICE_DEBUG("Map conservative");
}

void AxialGeoMultiscaleMapping::mapConsistent(const time::Sample &inData, Eigen::VectorXd &outData)
{
  PRECICE_TRACE();

  const int              inDataDimensions  = inData.dataDims;
  const Eigen::VectorXd &inputValues       = inData.values;
  Eigen::VectorXd       &outputValues      = outData;
  const int              outDataDimensions = outData.size() / output()->nVertices();

  PRECICE_ASSERT((inputValues.size() / static_cast<std::size_t>(inDataDimensions) == input()->nVertices()),
                 inputValues.size(), inDataDimensions, input()->nVertices());
  PRECICE_ASSERT((outputValues.size() / static_cast<std::size_t>(outDataDimensions) == output()->nVertices()),
                 outputValues.size(), outDataDimensions, output()->nVertices());

  // We currently don't support 1D data, so we need that the user specifies data of the same dimensions on both sides
  PRECICE_ASSERT(inDataDimensions == outDataDimensions);

  int effectiveCoordinate;
  if (inDataDimensions == 1) {
    effectiveCoordinate = 0;
  } else {
    effectiveCoordinate = static_cast<std::underlying_type_t<MultiscaleType>>(_axis);
    PRECICE_ASSERT(effectiveCoordinate == static_cast<std::underlying_type_t<MultiscaleType>>(MultiscaleAxis::X) ||
                       effectiveCoordinate == static_cast<std::underlying_type_t<MultiscaleType>>(MultiscaleAxis::Y) ||
                       effectiveCoordinate == static_cast<std::underlying_type_t<MultiscaleType>>(MultiscaleAxis::Z),
                   "Unknown multiscale axis type.");
  }

  PRECICE_DEBUG("Map consistent");

  if (_type == MultiscaleType::SPREAD) {
    size_t const inSize  = input()->nVertices();
    size_t const outSize = output()->nVertices();
    if (_dimension == MultiscaleDimension::D1D3 || _dimension == MultiscaleDimension::D1D2) {
      PRECICE_ASSERT(input()->nVertices() == 1);
      for (size_t i = 0; i < outSize; i++) {
        PRECICE_ASSERT(static_cast<size_t>((i * outDataDimensions) + effectiveCoordinate) < static_cast<size_t>(outputValues.size()), ((i * outDataDimensions) + effectiveCoordinate), outputValues.size());
        if (_profile == SpreadProfile::PARABOLIC) {
          if (_dimension == MultiscaleDimension::D1D2) {
            constexpr double factor                                     = 1.5;
            outputValues((i * outDataDimensions) + effectiveCoordinate) = factor * inputValues(effectiveCoordinate) * (1 - (_vertexDistances[i] * _vertexDistances[i]));
          } else {
            if (_crossSection == MultiscaleCrossSection::CIRCLE) {
              constexpr double factor                                     = 2.0;
              outputValues((i * outDataDimensions) + effectiveCoordinate) = factor * inputValues(effectiveCoordinate) * (1 - (_vertexDistances[i] * _vertexDistances[i]));
            } else if (_crossSection == MultiscaleCrossSection::SQUARE) {
              const double s1 = _vertexTransverseCoords[i][0];
              const double s2 = _vertexTransverseCoords[i][1];

              constexpr double factor = 2.096;
              constexpr double m      = 0.879;
              const double     b1raw  = 1.0 - s1 * s1;
              const double     b2raw  = 1.0 - s2 * s2;
              const double     b1     = std::max(0.0, b1raw);

              const double b2                                             = std::max(0.0, b2raw);
              outputValues((i * outDataDimensions) + effectiveCoordinate) = factor * inputValues(effectiveCoordinate) * std::pow(b1, m) * std::pow(b2, m);
            }
          }
        } else if (_profile == SpreadProfile::UNIFORM) {
          outputValues((i * outDataDimensions) + effectiveCoordinate) = inputValues(effectiveCoordinate);
        }
      }
    } else {
      PRECICE_ASSERT(_dimension == MultiscaleDimension::D2D3);
      PRECICE_ASSERT(_nearestVertex.size() == outSize);
      PRECICE_ASSERT(_vertexDistances.size() == outSize);
      PRECICE_ASSERT(_maxDistancePerInput.size() == inSize);
      for (size_t i = 0; i < outSize; i++) {
        PRECICE_ASSERT(_nearestVertex[i] >= 0 && static_cast<size_t>(_nearestVertex[i]) < inSize, _nearestVertex[i], inSize);
        PRECICE_ASSERT(static_cast<size_t>((i * outDataDimensions) + effectiveCoordinate) < static_cast<size_t>(outputValues.size()), (i * outDataDimensions) + effectiveCoordinate, outputValues.size());
        PRECICE_ASSERT(static_cast<size_t>((static_cast<size_t>(_nearestVertex[i]) * inDataDimensions) + effectiveCoordinate) < static_cast<size_t>(inputValues.size()), (static_cast<size_t>(_nearestVertex[i]) * inDataDimensions) + effectiveCoordinate, inputValues.size());
        double R = _maxDistancePerInput[static_cast<size_t>(_nearestVertex[i])];
        if (_profile == SpreadProfile::UNIFORM) {
          outputValues((i * outDataDimensions) + effectiveCoordinate) = inputValues((static_cast<size_t>(_nearestVertex[i]) * inDataDimensions) + effectiveCoordinate);
        } else if (_profile == SpreadProfile::PARABOLIC) {
          if (_crossSection == MultiscaleCrossSection::CIRCLE) {
            double r_hat                                                = _vertexDistances[i] / R;
            outputValues((i * outDataDimensions) + effectiveCoordinate) = (4.0 / 3.0) * inputValues((static_cast<size_t>(_nearestVertex[i]) * inDataDimensions) + effectiveCoordinate) * (1.0 - r_hat * r_hat);
          } else if (_crossSection == MultiscaleCrossSection::SQUARE) {
            constexpr double umax_over_umean = 2.096;
            constexpr double line_factor     = 1.5;
            constexpr double m               = 0.879;
            constexpr double eps             = 1e-12;
            const size_t     inIdx           = static_cast<size_t>(_nearestVertex[i]);
            const double     s1              = input()->vertex(inIdx).getCoords()[_lineCoord] / _radius;
            const double     s2              = _vertexDistances[i] / _radius;
            const double     b1raw           = 1.0 - s1 * s1;
            const double     b2raw           = 1.0 - s2 * s2;
            const double     b1              = std::max(eps, b1raw);

            const double b2                                             = std::max(0.0, b2raw);
            outputValues((i * outDataDimensions) + effectiveCoordinate) = inputValues((inIdx * inDataDimensions) + effectiveCoordinate) * (umax_over_umean / line_factor) * std::pow(b1, m - 1.0) * std::pow(b2, m);
          }
        }
      }
    }
  } else {
    PRECICE_ASSERT(_type == MultiscaleType::COLLECT);
    size_t const inSize  = input()->nVertices();
    size_t const outSize = output()->nVertices();
    if (_dimension == MultiscaleDimension::D1D3) {
      PRECICE_ASSERT(output()->nVertices() == 1);
      outputValues(effectiveCoordinate) = 0.0;
      for (size_t i = 0; i < inSize; i++) {
        PRECICE_ASSERT(static_cast<size_t>((i * inDataDimensions) + effectiveCoordinate) < static_cast<size_t>(inputValues.size()),
                       ((i * inDataDimensions) + effectiveCoordinate), inputValues.size());
        outputValues(effectiveCoordinate) +=
            inputValues((i * inDataDimensions) + effectiveCoordinate);
      }
      outputValues(effectiveCoordinate) = outputValues(effectiveCoordinate) / inSize;

    } else if (_dimension == MultiscaleDimension::D1D2) {
      PRECICE_ASSERT(output()->nVertices() == 1);

      outputValues(effectiveCoordinate) = 0.0;
      for (size_t i = 0; i < inSize; ++i) {
        PRECICE_ASSERT(static_cast<size_t>((i * inDataDimensions) + effectiveCoordinate) < static_cast<size_t>(inputValues.size()),
                       ((i * inDataDimensions) + effectiveCoordinate), inputValues.size());
        outputValues(effectiveCoordinate) +=
            inputValues((i * inDataDimensions) + effectiveCoordinate);
      }
      outputValues(effectiveCoordinate) = outputValues(effectiveCoordinate) / inSize;

    } else {
      PRECICE_ASSERT(_dimension == MultiscaleDimension::D2D3);
      PRECICE_ASSERT(_collectBands.size() == output()->nVertices());
      for (size_t j = 0; j < outSize; j++) {
        if (_collectBands[j].empty()) {
          continue;
        }
        PRECICE_ASSERT(static_cast<size_t>((j * outDataDimensions) + effectiveCoordinate) < static_cast<size_t>(outputValues.size()), (j * outDataDimensions) + effectiveCoordinate, outputValues.size());
        outputValues((j * outDataDimensions) + effectiveCoordinate) = 0.0;
        for (size_t k = 0; k < _collectBands[j].size(); k++) {
          int i = _collectBands[j][k];
          PRECICE_ASSERT(i >= 0 && static_cast<size_t>(i) < inSize, i, inSize);
          PRECICE_ASSERT(static_cast<size_t>((static_cast<size_t>(i) * inDataDimensions) + effectiveCoordinate) < static_cast<size_t>(inputValues.size()), (static_cast<size_t>(i) * inDataDimensions) + effectiveCoordinate, inputValues.size());
          outputValues((j * outDataDimensions) + effectiveCoordinate) += inputValues((static_cast<size_t>(i) * inDataDimensions) + effectiveCoordinate);
        }
        outputValues((j * outDataDimensions) + effectiveCoordinate) = outputValues((j * outDataDimensions) + effectiveCoordinate) / static_cast<double>(_collectBands[j].size());

        if (_profile == SpreadProfile::PARABOLIC) {
          if (_crossSection == MultiscaleCrossSection::CIRCLE) {
            outputValues((j * outDataDimensions) + effectiveCoordinate) *= 9.0 / 8.0;
          } else if (_crossSection == MultiscaleCrossSection::SQUARE) {
            const double s1    = output()->vertex(j).getCoords()[_lineCoord] / _radius;
            const double b1raw = 1.0 - s1 * s1;
            const double b1    = std::max(0.0, b1raw);
            outputValues((j * outDataDimensions) + effectiveCoordinate) *= 1.03745 * std::pow(b1, 0.121);
          }
        }
      }
    }
  }
}

void AxialGeoMultiscaleMapping::tagMeshFirstRound()
{
  PRECICE_TRACE();

  computeMapping();

  if (getConstraint() == CONSISTENT) {
    PRECICE_ASSERT(_type == MultiscaleType::SPREAD, "Not yet implemented");
    PRECICE_ASSERT(input()->nVertices() == 1);

    input()->tagAll();

  } else {
    PRECICE_ASSERT(getConstraint() == CONSERVATIVE);
    PRECICE_UNREACHABLE("Axial conservative geometric multiscale mapping is not implemented");
  }

  clear();
}

void AxialGeoMultiscaleMapping::tagMeshSecondRound()
{
  PRECICE_TRACE();
  // no operation needed here for the moment
}

std::string AxialGeoMultiscaleMapping::getName() const
{
  return "axial-geomultiscale";
}

} // namespace precice::mapping
