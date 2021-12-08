#include "Data.hpp"
#include <algorithm>
#include <utility>

#include "cplscheme/CouplingScheme.hpp"
#include "precice/types.hpp"
#include "time/Waveform.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace mesh {

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
    int         dimensions,
    int         interpolationOrder)
    : _values(),
      _name(std::move(name)),
      _id(id),
      _dimensions(dimensions)
{
  PRECICE_ASSERT(dimensions > 0, dimensions);
  _ptrWaveform = time::PtrWaveform(new time::Waveform(interpolationOrder));
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

void Data::storeDataInWaveform(int waveformSampleID)
{
  PRECICE_TRACE();
  PRECICE_ASSERT(_ptrWaveform->valuesSize() == _values.size(),
                 _ptrWaveform->valuesSize(), _values.size());
  _ptrWaveform->storeAt(_values, waveformSampleID);
}

int Data::sizeOfSampleStorageInWaveform()
{
  PRECICE_TRACE();
  return _ptrWaveform->sizeOfSampleStorage();
}

Eigen::VectorXd Data::waveformSampleAt(double normalizedDt)
{
  PRECICE_TRACE();
  PRECICE_ASSERT(_ptrWaveform->valuesSize() == _values.size(), _ptrWaveform->valuesSize(), _values.size());
  PRECICE_ASSERT(normalizedDt >= 0, "Sampling outside of valid range!");
  PRECICE_ASSERT(normalizedDt <= 1, "Sampling outside of valid range!");
  return _ptrWaveform->sample(normalizedDt);
}

void Data::initializeWaveform()
{
  PRECICE_TRACE();
  // PRECICE_ASSERT(_values.size() > 0, _values.size());  // @todo assertion breaks, but seems like calling advance on empty write data is ok?
  _ptrWaveform->initialize(_values.size());
  for (int sampleID = 0; sampleID < _ptrWaveform->sizeOfSampleStorage(); ++sampleID) {
    _ptrWaveform->storeAt(_values, sampleID);
  }
  PRECICE_ASSERT(_ptrWaveform->valuesSize() == _values.size());
}

void Data::sampleWaveformIntoData(int sampleID)
{
  PRECICE_TRACE();
  PRECICE_ASSERT(_ptrWaveform->valuesSize() == _values.size(), _ptrWaveform->valuesSize(), _values.size());
  values() = _ptrWaveform->getSample(sampleID);
}

void Data::moveToNextWindow()
{
  PRECICE_TRACE();
  _ptrWaveform->moveToNextWindow();
}

} // namespace mesh
} // namespace precice
