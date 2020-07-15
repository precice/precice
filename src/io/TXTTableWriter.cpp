#include "TXTTableWriter.hpp"
#include <algorithm>
#include <iomanip>
#include "logging/LogMacros.hpp"
#include "utils/Helpers.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace io {

TXTTableWriter::TXTTableWriter(
    const std::string &filename)
    : _data(),
      _writeIterator(_data.end()),
      _outputStream()
{
  _outputStream.open(filename);
  PRECICE_CHECK(_outputStream, "TXT table writer failed to open file \"" << filename << '"');

  _outputStream.setf(std::ios::showpoint);
  _outputStream.setf(std::ios::fixed);
  _outputStream << std::setprecision(16);
}

void TXTTableWriter::addData(
    const std::string &name,
    DataType           type)
{
  PRECICE_ASSERT(_outputStream);
  Data data;
  data.name = name;
  data.type = type;
  _data.push_back(data);
  if ((type == INT) || (type == DOUBLE)) {
    _outputStream << name << "  ";
  } else if (type == VECTOR2D) {
    for (int i = 0; i < 2; i++) {
      _outputStream << name << i << "  ";
    }
  } else {
    PRECICE_ASSERT(type == VECTOR3D);
    for (int i = 0; i < 3; i++) {
      _outputStream << name << i << "  ";
    }
  }
  _writeIterator = _data.end();
}

void TXTTableWriter::writeData(
    const std::string &name,
    int                value)
{
  PRECICE_ASSERT(_outputStream);
  PRECICE_ASSERT(not _data.empty());
  if (_writeIterator == _data.end()) {
    _writeIterator = _data.begin();
    _outputStream << "\n";
  }
  PRECICE_ASSERT(_writeIterator->name == name, _writeIterator->name, name);
  PRECICE_ASSERT(_writeIterator->type == INT, _writeIterator->type);
  _outputStream << value << "  ";
  _writeIterator++;
  if (_writeIterator == _data.end()) {
    _outputStream.flush();
  }
}

void TXTTableWriter::writeData(
    const std::string &name,
    double             value)
{
  PRECICE_ASSERT(_outputStream);
  PRECICE_ASSERT(not _data.empty());
  if (_writeIterator == _data.end()) {
    _writeIterator = _data.begin();
    _outputStream << "\n";
  }
  PRECICE_ASSERT(_writeIterator->name == name, _writeIterator->name, name);
  PRECICE_ASSERT(_writeIterator->type == DOUBLE, _writeIterator->type);
  _outputStream << value << "  ";
  _writeIterator++;
  if (_writeIterator == _data.end()) {
    _outputStream.flush();
  }
}

void TXTTableWriter::writeData(
    const std::string &    name,
    const Eigen::Vector2d &value)
{
  PRECICE_ASSERT(_outputStream);
  PRECICE_ASSERT(not _data.empty());
  if (_writeIterator == _data.end()) {
    _writeIterator = _data.begin();
    _outputStream << "\n";
  }
  PRECICE_ASSERT(_writeIterator->name == name, _writeIterator->name, name);
  PRECICE_ASSERT(_writeIterator->type == VECTOR2D, _writeIterator->type);
  for (int i = 0; i < value.size(); i++) {
    _outputStream << value[i] << "  ";
  }
  _writeIterator++;
  if (_writeIterator == _data.end()) {
    _outputStream.flush();
  }
}

void TXTTableWriter::writeData(
    const std::string &    name,
    const Eigen::Vector3d &value)
{
  PRECICE_ASSERT(_outputStream);
  PRECICE_ASSERT(not _data.empty());
  if (_writeIterator == _data.end()) {
    _writeIterator = _data.begin();
    _outputStream << "\n";
  }
  PRECICE_ASSERT(_writeIterator->name == name, _writeIterator->name, name);
  PRECICE_ASSERT(_writeIterator->type == VECTOR3D, _writeIterator->type);
  for (int i = 0; i < value.size(); i++) {
    _outputStream << value[i] << "  ";
  }
  _writeIterator++;
  if (_writeIterator == _data.end()) {
    _outputStream.flush();
  }
}

void TXTTableWriter::close()
{
  PRECICE_ASSERT(_outputStream.is_open());
  _outputStream.close();
}

/// Resets the table information.
void TXTTableWriter::reset()
{
  _data.clear();
  _writeIterator = _data.end();
}

} // namespace io
} // namespace precice
