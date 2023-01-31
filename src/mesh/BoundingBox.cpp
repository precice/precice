#include "BoundingBox.hpp"
#include <algorithm>
#include <limits>
#include <ostream>
#include <string>
#include <utility>
#include "logging/LogMacros.hpp"
#include "logging/Logger.hpp"
#include "mesh/Vertex.hpp"
#include "utils/assertion.hpp"

namespace precice::mesh {

logging::Logger BoundingBox::_log{"mesh::BoundingBox"};

BoundingBox::BoundingBox(Eigen::VectorXd boundMin, Eigen::VectorXd boundMax)
{
  PRECICE_ASSERT(boundMin.rows() == boundMax.rows(), "Dimension of min {} and max {} vertices should be the same.", boundMin.rows(), boundMax.rows());
  PRECICE_ASSERT((boundMin - boundMax).maxCoeff() < 0, "Each component of min vertex {} must be <= max vertex {} in the same axis direction.", boundMin, boundMax);

  _boundMin   = std::move(boundMin);
  _boundMax   = std::move(boundMax);
  _dimensions = _boundMin.rows();
}

BoundingBox::BoundingBox(std::vector<double> bounds)
{
  PRECICE_ASSERT((int) bounds.size() == 4 || (int) bounds.size() == 6, "Dimension of a bounding box can only be 2 or 3. The dimension is {}.", bounds.size() / 2);
  _dimensions = bounds.size() / 2;
  _boundMin   = Eigen::Map<Eigen::VectorXd, 0, Eigen::InnerStride<2>>(bounds.data(), _dimensions);
  _boundMax   = Eigen::Map<Eigen::VectorXd, 0, Eigen::InnerStride<2>>(bounds.data() + 1, _dimensions);
}

BoundingBox::BoundingBox(int dimension)
    : _dimensions(dimension)
{
  PRECICE_ASSERT(dimension == 2 || dimension == 3, "Dimension of a bounding box can only be 2 or 3 (in this case, it is {}).", dimension);

  _boundMin = Eigen::VectorXd(dimension);
  _boundMax = Eigen::VectorXd(dimension);

  // Define 'illegal' BoundingBox: _boundMin > _boundMax
  std::fill(_boundMin.data(), _boundMin.data() + _boundMin.rows(), std::numeric_limits<double>::max());
  std::fill(_boundMax.data(), _boundMax.data() + _boundMax.rows(), std::numeric_limits<double>::lowest());
}

bool BoundingBox::operator==(const BoundingBox &otherBB) const
{
  PRECICE_ASSERT(_dimensions == otherBB._dimensions, "Bounding boxes with different dimensions ({} and {}) cannot be compared.", _dimensions, otherBB._dimensions);

  return _boundMin.isApprox(otherBB._boundMin) && _boundMax.isApprox(otherBB._boundMax);
}

bool BoundingBox::empty() const
{
  return (_boundMax - _boundMin).isZero();
}

bool BoundingBox::isDefault() const
{
  return _boundMin.isApproxToConstant(std::numeric_limits<double>::max()) && _boundMax.isApproxToConstant(std::numeric_limits<double>::lowest());
}

bool BoundingBox::contains(const mesh::Vertex &vertex) const
{
  PRECICE_ASSERT(_dimensions == vertex.getDimensions(), "Vertex with different dimensions than this bounding box cannot be checked.");
  const auto &coords = vertex.rawCoords();
  for (int d = 0; d < _dimensions; d++) {
    if (coords[d] < _boundMin[d] || coords[d] > _boundMax[d]) {
      return false;
    }
  }
  return true;
}

Eigen::VectorXd BoundingBox::center() const
{
  PRECICE_ASSERT(!isDefault(), "Data of the bounding box is at default state.");
  Eigen::VectorXd cog(_dimensions);
  for (int d = 0; d < _dimensions; d++) {
    cog[d] = (_boundMax[d] + _boundMin[d]) / 2.0;
  }
  return cog;
}

Eigen::VectorXd BoundingBox::minCorner() const
{
  return _boundMin;
}

Eigen::VectorXd BoundingBox::maxCorner() const
{
  return _boundMax;
}

double BoundingBox::getArea(std::vector<bool> deadAxis)
{
  PRECICE_ASSERT(!isDefault(), "Data of the bounding box is at default state.");
  double meshArea = 1.0;
  for (int d = 0; d < _dimensions; d++)
    if (not deadAxis[d])
      meshArea *= _boundMax[d] - _boundMin[d];
  return meshArea;
}

int BoundingBox::getDimension() const
{
  return _dimensions;
}

const std::vector<double> BoundingBox::dataVector() const
{
  std::vector<double> bounds({_boundMin[0], _boundMax[0], _boundMin[1], _boundMax[1]});
  if (_dimensions == 3) {
    bounds.insert(bounds.end(), {_boundMin[2], _boundMax[2]});
  }
  return bounds;
}

void BoundingBox::expandBy(const BoundingBox &otherBB)
{
  PRECICE_ASSERT(_dimensions == otherBB.getDimension(), "Other BoundingBox with different dimensions than this bounding box cannot be used to expand it.");
  for (int d = 0; d < _dimensions; d++) {
    _boundMin[d] = std::min(_boundMin[d], otherBB._boundMin[d]);
    _boundMax[d] = std::max(_boundMax[d], otherBB._boundMax[d]);
  }
}

void BoundingBox::expandBy(const Vertex &vertices)
{
  PRECICE_ASSERT(_dimensions == vertices.getDimensions(), "Vertex with different dimensions than this bounding box cannot be used to expand it.");
  const auto coords = vertices.rawCoords();
  for (int d = 0; d < _dimensions; ++d) {
    _boundMin[d] = std::min(coords[d], _boundMin[d]);
    _boundMax[d] = std::max(coords[d], _boundMax[d]);
  }
}

void BoundingBox::expandBy(double value)
{
  if (!isDefault()) {
    for (int d = 0; d < _dimensions; d++) {
      _boundMin[d] -= value;
      _boundMax[d] += value;
    }
  }
}

void BoundingBox::scaleBy(double safetyFactor)
{
  if (!isDefault()) {
    double maxSideLength = 1e-6; // we need some minimum > 0 here
    maxSideLength        = (_boundMax - _boundMin).maxCoeff();
    _boundMax.array() += safetyFactor * maxSideLength;
    _boundMin.array() -= safetyFactor * maxSideLength;
    PRECICE_DEBUG("Merged BoundingBox {}", *this);
  }
}

bool BoundingBox::overlapping(const BoundingBox &otherBB)
{
  for (int d = 0; d < _dimensions; d++) {
    if ((_boundMin[d] < otherBB._boundMin[d] && _boundMax[d] < otherBB._boundMin[d]) ||
        (otherBB._boundMin[d] < _boundMin[d] && otherBB._boundMax[d] < _boundMin[d])) {
      return false;
    }
  }
  return true;
}

void BoundingBox::print(std::ostream &out) const
{
  out << "( ";
  for (int d = 0; d < _dimensions; ++d) {
    out << "[" << _boundMin[d] << " " << _boundMax[d] << "], ";
  }
  out << ")";
}

std::ostream &operator<<(std::ostream &os, const BoundingBox &bb)
{
  bb.print(os);
  return os;
}

} // namespace precice::mesh
