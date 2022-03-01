#include "Data.hpp"
#include <algorithm>
#include <utility>

#include "SharedPointer.hpp"
#include "precice/types.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace mesh {

Data::Data()
    : _name(""),
      _id(-1),
      _dimensions(0),
      _spacialDimensions(-1),
      _hasGradient(false)
{
  PRECICE_ASSERT(false);
}

Data::Data(
    std::string name,
    DataID      id,
    int         dimensions,
    int         spacialDimensions,
    bool        hasGradient)
    : _values(),
      _name(std::move(name)),
      _id(id),
      _dimensions(dimensions),
      _spacialDimensions(spacialDimensions),
      _hasGradient(hasGradient)
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
  auto begin = _values.data();
  auto end   = begin + _values.size();
  std::fill(begin, end, 0.0);

  if (_hasGradient) {
    auto beginGradient = _gradientValues.data();
    auto endGradient   = beginGradient + _gradientValues.size();
    std::fill(beginGradient, endGradient, 0.0);
  }
}

bool Data::hasGradient() const
{
  return _hasGradient;
}

int Data::getDimensions() const
{
  return _dimensions;
}

int Data::getSpacialDimensions() const
{
  return _spacialDimensions;
}

} // namespace mesh
} // namespace precice
