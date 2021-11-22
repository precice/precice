#include "Data.hpp"
#include <algorithm>
#include <utility>

#include "cplscheme/CouplingScheme.hpp"
#include "precice/types.hpp"
#include "time/Waveform.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace mesh {

// const int Data::EXTRAPOLATION_ORDER = cplscheme::CouplingScheme::UNDEFINED_EXTRAPOLATION_ORDER;
const int Data::EXTRAPOLATION_ORDER = 0; // @todo should be cplscheme::CouplingScheme::UNDEFINED_EXTRAPOLATION_ORDER, but then some tests break.
const int Data::INTERPOLATION_ORDER = 1;

size_t Data::_dataCount = 0;

Data::Data()
    : _name(""),
      _id(-1),
      _dimensions(0)
{
  PRECICE_ASSERT(false);
}

Data::Data(
    std::string name,
    DataID      id,
    int         dimensions)
    : _values(),
      _name(std::move(name)),
      _id(id),
      _dimensions(dimensions)
{
  PRECICE_ASSERT(dimensions > 0, dimensions);
  _ptrWaveform = time::PtrWaveform(new time::Waveform(Data::EXTRAPOLATION_ORDER, Data::INTERPOLATION_ORDER));
  _dataCount++;
}

Data::~Data()
{
  _dataCount--;
}

Eigen::VectorXd &Data::values()
{
  return _values;
}

const Eigen::VectorXd &Data::values() const
{
  return _values;
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
}

int Data::getDimensions() const
{
  return _dimensions;
}

size_t Data::getDataCount()
{
  return _dataCount;
}

void Data::resetDataCount()
{
  _dataCount = 0;
}

time::PtrWaveform Data::waveform()
{
  return _ptrWaveform;
}

void Data::setExtrapolationOrder(int extrapolationOrder)
{
  _ptrWaveform->setExtrapolationOrder(extrapolationOrder);
}

} // namespace mesh
} // namespace precice
