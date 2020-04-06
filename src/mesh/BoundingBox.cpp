#include "BoundingBox.hpp"

namespace precice {
namespace mesh {

BoundingBox::BoundingBox(std::vector<double> bounds, double safetyFactor)
{

  PRECICE_CHECK(bounds.size() == 4 || bounds.size() == 6, "Dimension of the bounding box should be 2 or 3.");
  _dimensions = bounds.size() / 2;

  for (int i = 0; i < _dimensions; ++i) {
    _bounds.at(i) = bounds.at(i);
  }
}

BoundingBox::BoundingBox(const BoundingBox &bb)
{
  _dimensions = bb._dimensions;
  _bounds     = bb._bounds;
}

BoundingBox::~BoundingBox()
{
}

std::ostream &operator<<(std::ostream &out, const BoundingBox &bb)
{
  for (int d = 0; d < bb._dimensions; ++d) {
    out << "dim: " << d << " min: " << bb._bounds[2 * d] << ", max: " << bb._bounds[2 * d + 1] << "\n";
  }
  return out;
}

BoundingBox BoundingBox::createFromData(std::vector<double> bounds, double safetyFactor)
{
  BoundingBox box{bounds, safetyFactor};
  return box;
}

void BoundingBox::setBounds(int dimension, double min, double max)
{
  if (_dimensions == 2) {
    PRECICE_CHECK(dimension == 0 || dimension == 1, "Given bound limit axis is not compatible with bounding box dimensions.");
  }
  if (_dimensions == 3) {
    PRECICE_CHECK(dimension == 0 || dimension == 1 || dimension == 2, "Given bound limit axis is not compatible with bounding box dimensions.");
  }

  _bounds.at(dimension * 2)     = min;
  _bounds.at(dimension * 2 + 1) = max;
}

bool BoundingBox::isVertexInBB(const mesh::Vertex &vertex)
{
  for (int d = 0; d < _dimensions; d++) {
    if (vertex.getCoords()[d] < _bounds[2 * d] || vertex.getCoords()[d] > _bounds[2 * d + 1]) {
      return false;
    }
  }
  return true;
}

void BoundingBox::mergeBoundingBoxes(const BoundingBox &otherBB, double safetyFactor)
{
  if (_prepared)
    return;

  PRECICE_ASSERT(!otherBB._bounds.empty(), "Output Mesh of from Mapping has an empty bounding box!");
  for (int d = 0; d < _dimensions; d++) {
    _bounds[2 * d]     = std::min(_bounds[2 * d], otherBB._bounds[2 * d]);
    _bounds[2 * d + 1] = std::max(_bounds[2 * d + 1], otherBB._bounds[2 * d + 1]);
  }

  // Enlarge BB
  PRECICE_ASSERT(safetyFactor >= 0.0);

  double maxSideLength = 1e-6; // we need some minimum > 0 here

  for (int d = 0; d < _dimensions; d++) {
    if (_bounds[2 * d + 1] > _bounds[2 * d])
      maxSideLength = std::max(maxSideLength, _bounds[2 * d + 1] - _bounds[2 * d]);
  }
  for (int d = 0; d < _dimensions; d++) {
    _bounds[2 * d + 1] += safetyFactor * maxSideLength;
    _bounds[2 * d] -= safetyFactor * maxSideLength;
    PRECICE_DEBUG("Merged BoundingBox, dim: " << d << ", first: " << _bounds[2 * d] << ", second: " << _bounds[2 * d + 1]);
  }

  _prepared = true;
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

} // namespace mesh
} // namespace precice