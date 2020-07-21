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

namespace precice {
namespace mesh {

logging::Logger BoundingBox::_log{"mesh::BoundingBox"};

BoundingBox::BoundingBox(std::vector<double> bounds)
{
  PRECICE_ASSERT((int) bounds.size() == 4 || (int) bounds.size() == 6, "Dimension of a bounding box can only be 2 or 3. Given dimension is: " << bounds.size() / 2);
  _bounds     = std::move(bounds);
  _dimensions = _bounds.size() / 2;
}

BoundingBox::BoundingBox(int dimension)
    : _dimensions(dimension)
{
  PRECICE_ASSERT(dimension == 2 || dimension == 3, "Dimension of a bounding box can only be 2 or 3. Given dimension is: " << dimension);
  for (int i = 0; i < _dimensions; ++i) {
    _bounds.push_back(std::numeric_limits<double>::max());
    _bounds.push_back(std::numeric_limits<double>::lowest());
  }
}

bool BoundingBox::operator==(const BoundingBox &otherBB) const
{
  PRECICE_ASSERT(_dimensions == otherBB._dimensions, "Bounding boxes with different dimensions cannot be compared. Dimensions: " << _dimensions << " and " << otherBB._dimensions);
  for (int i = 0; i < _dimensions; ++i) {
    if (_bounds.at(i) != otherBB._bounds.at(i)) {
      return false;
    }
  }
  return true;
}

bool BoundingBox::empty() const
{
  for (int i = 0; i < _dimensions; ++i) {
    if (_bounds[2 * i] > _bounds[2 * i + 1]) {
      return true;
    }
  }
  return false;
}

bool BoundingBox::contains(const mesh::Vertex &vertex) const
{
  PRECICE_ASSERT(_dimensions == vertex.getDimensions(), "Vertex with different dimensions than bounding box cannot be checked.");
  for (int d = 0; d < _dimensions; d++) {
    if (vertex.getCoords()[d] < _bounds.at(2 * d) || vertex.getCoords()[d] > _bounds.at(2 * d + 1)) {
      return false;
    }
  }
  return true;
}

Eigen::VectorXd BoundingBox::center() const
{
  PRECICE_ASSERT(!empty(), "Data of the bounding box is at default state.");
  Eigen::VectorXd cog(_dimensions);
  for (int d = 0; d < _dimensions; d++) {
    cog[d] = (_bounds[2 * d + 1] - _bounds[2 * d]) / 2.0 + _bounds[2 * d];
  }
  return cog;
}

Eigen::VectorXd BoundingBox::minCorner() const
{
  Eigen::VectorXd min(_dimensions);
  for (int d = 0; d < _dimensions; d++) {
    min[d] = _bounds[2 * d];
  }
  return min;
}

Eigen::VectorXd BoundingBox::maxCorner() const
{
  Eigen::VectorXd max(_dimensions);
  for (int d = 0; d < _dimensions; d++) {
    max[d] = _bounds[2 * d + 1];
  }
  return max;
}

double BoundingBox::getArea(std::vector<bool> deadAxis)
{
  PRECICE_ASSERT(!empty(), "Data of the bounding box is at default state.");
  double meshArea = 1.0;
  for (int d = 0; d < _dimensions; d++)
    if (not deadAxis[d])
      meshArea *= _bounds[2 * d + 1] - _bounds[2 * d];
  return meshArea;
}

int BoundingBox::getDimension() const
{
  return _dimensions;
}

const std::vector<double> BoundingBox::dataVector() const
{
  return _bounds;
}

void BoundingBox::expandBy(const BoundingBox &otherBB)
{
  for (int d = 0; d < _dimensions; d++) {
    _bounds[2 * d]     = std::min(_bounds[2 * d], otherBB._bounds[2 * d]);
    _bounds[2 * d + 1] = std::max(_bounds[2 * d + 1], otherBB._bounds[2 * d + 1]);
  }
}

void BoundingBox::expandBy(const Vertex &vertices)
{
  PRECICE_ASSERT(_dimensions == vertices.getDimensions(), "Vertex with different dimensions than bounding box cannot be used to expand bounding box");
  for (int d = 0; d < _dimensions; ++d) {
    _bounds.at(2 * d)     = std::min(vertices.getCoords()[d], _bounds.at(2 * d));
    _bounds.at(2 * d + 1) = std::max(vertices.getCoords()[d], _bounds.at(2 * d + 1));
  }
}

void BoundingBox::expandBy(double value)
{
  if (empty())
    return;
  for (int d = 0; d < _dimensions; d++) {
    _bounds[2 * d] -= value;
    _bounds[2 * d + 1] += value;
  }
}

void BoundingBox::scaleBy(double safetyFactor)
{
  if (empty())
    return;
  double maxSideLength = 1e-6; // we need some minimum > 0 here
  for (int d = 0; d < _dimensions; d++) {
    if (_bounds.at(2 * d + 1) > _bounds.at(2 * d))
      maxSideLength = std::max(maxSideLength, _bounds[2 * d + 1] - _bounds[2 * d]);
  }
  for (int d = 0; d < _dimensions; d++) {
    _bounds.at(2 * d + 1) += safetyFactor * maxSideLength;
    _bounds.at(2 * d) -= safetyFactor * maxSideLength;
    PRECICE_DEBUG("Merged BoundingBox" << *this);
  }
}

bool BoundingBox::overlapping(const BoundingBox &otherBB)
{
  for (int d = 0; d < _dimensions; d++) {
    if ((_bounds[2 * d] < otherBB._bounds[2 * d] && _bounds[2 * d + 1] < otherBB._bounds[2 * d]) ||
        (otherBB._bounds[2 * d] < _bounds[2 * d] && otherBB._bounds[2 * d + 1] < _bounds[2 * d])) {
      return false;
    }
  }
  return true;
}

void BoundingBox::print(std::ostream &out) const
{
  out << "( ";
  for (int d = 0; d < _dimensions; ++d) {
    out << "[" << _bounds[2 * d] << " " << _bounds[2 * d + 1] << "], ";
  }
  out << ")";
}

std::ostream &operator<<(std::ostream &os, const BoundingBox &bb)
{
  bb.print(os);
  return os;
}

} // namespace mesh
} // namespace precice
