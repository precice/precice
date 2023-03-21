#include "GlobalData.hpp"
#include <algorithm>
#include <utility>

#include "SharedPointer.hpp"
#include "precice/types.hpp"
#include "utils/assertion.hpp"

namespace precice::mesh {

GlobalData::GlobalData()
    : _name(""),
      _id(-1),
      _dimensions(0),
      _spatialDimensions(-1)
{
  PRECICE_ASSERT(false);
}

GlobalData::GlobalData(
    std::string name,
    DataID      id,
    int         dimensions,
    int         spatialDimensions)
    : _values(dimensions),
      _name(std::move(name)),
      _id(id),
      _dimensions(dimensions),
      _spatialDimensions(spatialDimensions)
{
  PRECICE_ASSERT(dimensions > 0, dimensions);
}

Eigen::VectorXd &GlobalData::values()
{
  return _values;
}

const Eigen::VectorXd &GlobalData::values() const
{
  return _values;
}

const std::string &GlobalData::getName() const
{
  return _name;
}

DataID GlobalData::getID() const
{
  return _id;
}

void GlobalData::toZero()
{
  _values.setZero();
}

int GlobalData::getDimensions() const
{
  return _dimensions;
}

int GlobalData::getSpatialDimensions() const
{
  return _spatialDimensions;
}

} // namespace precice::mesh
