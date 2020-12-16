#include "Gradient.hpp"
#include <algorithm>
#include "utils/assertion.hpp"

namespace precice {
namespace mesh {

size_t Gradient::_gradientCount = 0;

Gradient::Gradient()
    : _name(""),
      _id(-1),
      _dimensions(0)
{
  PRECICE_ASSERT(false);
}

Gradient::Gradient(
    const std::string &name,
    int                id,
    int                dimensions)
    : _values(),
      _name(name),
      _id(id),
      _dimensions(dimensions)
{
  PRECICE_ASSERT(dimensions > 0, dimensions);
  _values.resize(dimensions,0);
  PRECICE_ASSERT(_values.rows() == dimensions, _values.rows(), dimensions);
  _gradientCount++;
}

Gradient::~Gradient()
{
  _gradientCount--;
}

Eigen::MatrixXd &Gradient::values()
{
  return _values;
}

const Eigen::MatrixXd &Gradient::values() const
{
  return _values;
}

const std::string &Gradient::getName() const
{
  return _name;
}

int Gradient::getID() const
{
  return _id;
}

void Gradient::toZero()
{
  auto begin = _values.data();
  auto end   = begin + _values.size();
  std::fill(begin, end, 0.0);
}

int Gradient::getDimensions() const
{
  return _dimensions;
}

size_t Gradient::getGradientCount()
{
  return _gradientCount;
}

void Gradient::resetGradientCount()
{
  _gradientCount = 0;
}

} // namespace mesh
} // namespace precice
