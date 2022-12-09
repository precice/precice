#include "Data.hpp"
#include <algorithm>
#include <utility>

#include "SharedPointer.hpp"
#include "precice/types.hpp"
#include "utils/assertion.hpp"

namespace precice::mesh {

Data::Data()
    : _name(""),
      _id(-1),
      _dimensions(0),
      _spatialDimensions(-1)
{
  PRECICE_ASSERT(false);
}

Data::Data(
    std::string name,
    DataID      id,
    int         dimensions,
    int         spatialDimensions)
    : _values(),
      _name(std::move(name)),
      _id(id),
      _dimensions(dimensions),
      _spatialDimensions(spatialDimensions)
{
  PRECICE_ASSERT(dimensions > 0, dimensions);
}

Eigen::VectorXd &Data::values()
{
  return _values;
}

const Eigen::VectorXd &Data::values() const
{
  return _values;
}

Eigen::MatrixXd &Data::gradientValues()
{
  return _gradientValues;
}

const Eigen::MatrixXd &Data::gradientValues() const
{
  return _gradientValues;
}

const std::string &Data::getName() const
{
  return _name;
}

DataID Data::getID() const
{
  return _id;
}

void Data::toZero()
{
  _values.setZero();
  if (_hasGradient) {
    _gradientValues.setZero();
  }
}

bool Data::hasGradient() const
{
  return _hasGradient;
}

void Data::requireDataGradient()
{
  _hasGradient = true;
};

int Data::getDimensions() const
{
  return _dimensions;
}

int Data::getSpatialDimensions() const
{
  return _spatialDimensions;
}

} // namespace precice::mesh
