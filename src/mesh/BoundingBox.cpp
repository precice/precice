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
  // PRECICE_ASSERT 1
  const int &dimMin = boundMin.rows();
  const int &dimMax = boundMax.rows();
  PRECICE_ASSERT(dimMin == dimMax, "Dimension of min and max vertices should be the same.", dimMin, dimMax);

  // PRECICE_ASSERT 2
  bool isMinLessThanMax = true;
  for (int i = 0; i < dimMin; i++) {
    if (boundMin[i] >= boundMax[i]) {
      isMinLessThanMax = false;
      break;
    }
  }
  PRECICE_ASSERT(isMinLessThanMax, "Each component of min vertex must be less than that of max vertex in the same axis direction.", boundMin, boundMax);

  // Assign private members
  _boundMin   = std::move(boundMin); // converts lval to rval
  _boundMax   = std::move(boundMax);
  _dimensions = dimMin;
}

BoundingBox::BoundingBox(std::vector<double> bounds)
{
  PRECICE_ASSERT((int) bounds.size() == 4 || (int) bounds.size() == 6, "Dimension of a bounding box can only be 2 or 3.", bounds.size() / 2);
  _dimensions = _bounds.size() / 2;
  _boundMin   = Eigen::Map<Eigen::VectorXd>(bounds.data(), _dimensions);
  _boundMax   = Eigen::Map<Eigen::VectorXd>(bounds.data() + _dimensions, _dimensions);
}

// TODO: use more proper way to assign _boundMin & _boundMax
BoundingBox::BoundingBox(int dimension)
    : _dimensions(dimension)
{
  PRECICE_ASSERT(dimension == 2 || dimension == 3, "Dimension of a bounding box can only be 2 or 3.", dimension);
  double min = std::numeric_limits<double>::lowest();
  double max = std::numeric_limits<double>::max();

  if (dimension == 2) {
    _boundMin = Eigen::VectorXd{min,min};
    _boundMax << max, max;
  }
  else {
    _boundMin << min, min, min;
    _boundMax << max, max, max;
  }
}

bool BoundingBox::operator==(const BoundingBox &otherBB) const
{
  PRECICE_ASSERT(_dimensions == otherBB._dimensions, "Bounding boxes with different dimensions cannot be compared.", _dimensions, otherBB._dimensions);
  for (int i = 0; i < _dimensions; ++i) {
    if (_boundMin[i] != otherBB._boundMin[i] && _boundMax[i] != otherBB._boundMax[i]) {
      return false;
    }
  }
  return true;
}

bool BoundingBox::empty() const
{
  for (int i = 0; i < _dimensions; ++i) {
    if (_boundMin[i] > _boundMax[i]) {
      return true;
    }
  }
  return false;
}

bool BoundingBox::contains(const mesh::Vertex &vertex) const
{
  PRECICE_ASSERT(_dimensions == vertex.getDimensions(), "Vertex with different dimensions than bounding box cannot be checked.");
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
  PRECICE_ASSERT(!empty(), "Data of the bounding box is at default state.");
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
  PRECICE_ASSERT(!empty(), "Data of the bounding box is at default state.");
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
  std::vector<double> _bounds{_boundMin[0], _boundMax[0], _boundMin[1], _boundMax[1], _boundMin[2], _boundMax[2]};
  return _bounds;
}

void BoundingBox::expandBy(const BoundingBox &otherBB)
{
  for (int d = 0; d < _dimensions; d++) {
    _boundMin[d] = std::min(_boundMin[d], otherBB._boundMin[d]);
    _boundMax[d] = std::max(_boundMax[d], otherBB._boundMax[d]);
  }
}

void BoundingBox::expandBy(const Vertex &vertices)
{
  PRECICE_ASSERT(_dimensions == vertices.getDimensions(), "Vertex with different dimensions than bounding box cannot be used to expand bounding box");
  const auto coords = vertices.rawCoords();
  for (int d = 0; d < _dimensions; ++d) {
    _boundMin[d] = std::min(coords[d], _boundMin[d]);
    _boundMax[d] = std::max(coords[d], _boundMax[d]);
  }
}
// TODO: empty() definition
void BoundingBox::expandBy(double value)
{
  if (empty())
    return;
  for (int d = 0; d < _dimensions; d++) {
    _boundMin[d] -= value;
    _boundMax[d] += value;
  }
}
// TODO: empty() definition
void BoundingBox::scaleBy(double safetyFactor)
{
  if (empty())
    return;
  double maxSideLength = 1e-6; // we need some minimum > 0 here
  for (int d = 0; d < _dimensions; d++) {
    if (_boundMax[d] > _boundMin[d])
      maxSideLength = std::max(maxSideLength, _boundMax[d] - _boundMin[d]);
  }
  for (int d = 0; d < _dimensions; d++) {
    _boundMax[d] += safetyFactor * maxSideLength;
    _boundMin[d] -= safetyFactor * maxSideLength;
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
