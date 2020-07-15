#include "TXTReader.hpp"
#include "logging/LogMacros.hpp"

namespace precice {
namespace io {

TXTReader::TXTReader(
    const std::string &filename)
    : _file()
{
  _file.open(filename);
  PRECICE_CHECK(_file, "TXT reader failed to open file \"" << filename << '"');
  _file.setf(std::ios::showpoint);
  _file.setf(std::ios::fixed);
  //_file << std::setprecision(16);
}

} // namespace io
} // namespace precice
