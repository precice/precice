#include "BoundingBox.hpp"

namespace precice {
namespace mesh {

BoundingBox::BoundingBox(std::vector<double> bounds)
{
  //PRECICE_CHECK(bounds.size() == 4 || bounds.size() == 6, "Dimension of the bounding box should be 2 or 3.");
  _bounds = bounds;
  _dimensions = _bounds.size() / 2;
}

BoundingBox::BoundingBox(int dimension)
:_dimensions(dimension)
{
  PRECICE_CHECK(dimension == 2 || dimension == 3, "Dimension of the bounding box should be 2 or 3.");
  for(int i = 0; i < _dimensions; ++i){
    _bounds.push_back(std::numeric_limits<double>::max());
    _bounds.push_back(std::numeric_limits<double>::lowest());
  }
}

BoundingBox::BoundingBox(){
  _dimensions = 3;
  for(int i = 0; i < _dimensions; ++i){
    _bounds.push_back(std::numeric_limits<double>::max());
    _bounds.push_back(std::numeric_limits<double>::lowest());
  }
}

void BoundingBox::modifyForTest(int rank, std::string testName){
  if(testName == "com")
  {
    for (int i = 0; i < _dimensions; i++) {
      _bounds[2*i] = rank*i;
      _bounds[2*i + 1] = i + 1;
    }
  }
  if(testName == "partition"){
    for (int i = 0; i < _dimensions; i++) {
      _bounds[2*i] = 3 - rank - 1;
      _bounds[2*i + 1] = 3 - rank;
    }
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

bool BoundingBox::operator==(const BoundingBox& otherBB) const{
  if(_dimensions != otherBB._dimensions){
    return false;
  }
  for(int i = 0; i < _dimensions; ++i){
    if(std::abs(_bounds.at(i) - otherBB._bounds.at(i)) > 1e-6){
      return false;
    }
  }
  return true;
}

BoundingBox BoundingBox::createFromData(std::vector<double> bounds)
{
  BoundingBox box{bounds};
  return box;
}

void BoundingBox::setMin(int dimension, double min){
  PRECICE_CHECK(dimension == 2 || dimension == 3, "Given dimension to set bounding box is not compatible with the bounding box.");
  _bounds.at(2*dimension) = min;
}

void BoundingBox::setMax(int dimension, double max){
  PRECICE_CHECK(dimension == 2 || dimension == 3, "Given dimension to set bounding box is not compatible with the bounding box.");
  _bounds.at(2*dimension+1) = max;
}

void BoundingBox::setSafetyFactor(double safetyFactor){
  _safetyFactor = safetyFactor;
}

bool BoundingBox::isVertexInBB(const mesh::Vertex &vertex)
{
  for (int d = 0; d < _dimensions; d++) {
    if (vertex.getCoords()[d] < _bounds.at(2 * d) || vertex.getCoords()[d] > _bounds.at(2 * d + 1)) {
      return false;
    }
  }
  return true;
}

double BoundingBox::getData(int dimension, int type) const{
  if(type == 1){
    return _bounds.at(2*dimension);
  }
  if(type == 2){
    return _bounds.at(2*dimension+1);
  }
}

int BoundingBox::getDimension() const{
  return _dimensions;
}

int BoundingBox::getSize() const{
  return _bounds.size();
}

const double* BoundingBox::data() const{
  return _bounds.data();
}

const std::vector<double> BoundingBox::dataVector() const{
  return _bounds;
}

bool BoundingBox::mergeBoundingBoxes(const BoundingBox &otherBB)
{
  for (int d = 0; d < _dimensions; d++) {
    _bounds[2 * d]     = std::min(_bounds[2 * d], otherBB._bounds[2 * d]);
    _bounds[2 * d + 1] = std::max(_bounds[2 * d + 1], otherBB._bounds[2 * d + 1]);
  }

  // Enlarge BB
  PRECICE_ASSERT(_safetyFactor >= 0.0);

  double maxSideLength = 1e-6; // we need some minimum > 0 here

  for (int d = 0; d < _dimensions; d++) {
    if (_bounds.at(2 * d + 1) > _bounds.at(2 * d))
      maxSideLength = std::max(maxSideLength, _bounds[2 * d + 1] - _bounds[2 * d]);
  }
  for (int d = 0; d < _dimensions; d++) {
    _bounds.at(2 * d + 1) += _safetyFactor * maxSideLength;
    _bounds.at(2 * d) -= _safetyFactor * maxSideLength;
    PRECICE_DEBUG("Merged BoundingBox" << *this);
  }
  return true;
}

void BoundingBox::expandTo(const Vertex& vertices){

  for(int d = 0; d < _dimensions; ++d){
    _bounds.at(2*d) = std::min(vertices.getCoords()[d], _bounds.at(2*d));
    _bounds.at(2*d+1) = std::max(vertices.getCoords()[d], _bounds.at(2*d+1));
  }

}

void BoundingBox::enlargeWith(double value){
  for (int d = 0; d < _dimensions; d++) {
    _bounds[2*d] -= value;
    _bounds[2*d + 1] += value;   
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