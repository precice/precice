#include "BoundingBox.hpp"

namespace precice {
namespace mesh {

BoundingBox::BoundingBox(std::vector<double> bounds)
: _bounds(bounds), _dimensions(bounds.size()/2)
{

}

BoundingBox::BoundingBox(int dimension)
:_dimensions(dimension)
{
  for(int i = 0; i < _dimensions; ++i){
    _bounds.push_back(std::numeric_limits<double>::max());
    _bounds.push_back(std::numeric_limits<double>::lowest());
  }
}

BoundingBox::BoundingBox(const BoundingBox &bb)
{
  _dimensions = bb._dimensions;
  _bounds     = bb._bounds;
}

BoundingBox::BoundingBox(){}

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

bool BoundingBox::operator==(const BoundingBox& otherBB) const{
  return _bounds == otherBB._bounds;
}

BoundingBox BoundingBox::createFromData(std::vector<double> bounds)
{
  BoundingBox box{bounds};
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

void BoundingBox::setMin(int dimension, double min){
  _bounds[2*dimension] = min;
}

void BoundingBox::setMax(int dimension, double max){
  _bounds[2*dimension+1] = max;
}

void BoundingBox::setSafetyFactor(double safetyFactor){
  _safetyFactor = safetyFactor;
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

double BoundingBox::getData(int dimension, int type) const{

  if(type == 1){
    return _bounds[2*dimension];
  }
  if(type == 2){
    return _bounds[2*dimension+1];
  }

}

int BoundingBox::getDimension() const{
  return _dimensions;
}

const double* BoundingBox::data() const{

  return _bounds.data();

}

std::vector<double> BoundingBox::dataVector(){
  return _bounds;
}

bool BoundingBox::empty(){
  return _bounds.empty();
}

bool BoundingBox::mergeBoundingBoxes(const BoundingBox &otherBB)
{
  if (_prepared)
    return _prepared;

  PRECICE_ASSERT(!otherBB._bounds.empty(), "Output Mesh of from Mapping has an empty bounding box!");
  for (int d = 0; d < _dimensions; d++) {
    _bounds[2 * d]     = std::min(_bounds[2 * d], otherBB._bounds[2 * d]);
    _bounds[2 * d + 1] = std::max(_bounds[2 * d + 1], otherBB._bounds[2 * d + 1]);
  }

  // Enlarge BB
  PRECICE_ASSERT(_safetyFactor >= 0.0);

  double maxSideLength = 1e-6; // we need some minimum > 0 here

  for (int d = 0; d < _dimensions; d++) {
    if (_bounds[2 * d + 1] > _bounds[2 * d])
      maxSideLength = std::max(maxSideLength, _bounds[2 * d + 1] - _bounds[2 * d]);
  }
  for (int d = 0; d < _dimensions; d++) {
    _bounds[2 * d + 1] += _safetyFactor * maxSideLength;
    _bounds[2 * d] -= _safetyFactor * maxSideLength;
    PRECICE_DEBUG("Merged BoundingBox, dim: " << d << ", first: " << _bounds[2 * d] << ", second: " << _bounds[2 * d + 1]);
  }

  _prepared = true;
  return _prepared;
}

void BoundingBox::expandTo(const Vertex& vertices){

  for(int d = 0; d < _dimensions; ++d){
    _bounds[2*d] = std::min(vertices.getCoords()[d], _bounds[2*d]);
    _bounds[2*d+1] = std::max(vertices.getCoords()[d], _bounds[2*d+1]);
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

} // namespace mesh
} // namespace precice