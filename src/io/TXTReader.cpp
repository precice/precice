#include "TXTReader.hpp"
#include "logging/Logger.hpp"

namespace precice {
namespace io {

TXTReader:: TXTReader
(
  const std::string& filename )
:
  _file()
{
  _file.open(filename.c_str());
  if (not _file){
    ERROR("Could not open file \"" << filename << "\" for txt reading!");
  }
  _file.setf(std::ios::showpoint);
  _file.setf(std::ios::fixed);
  //_file << std::setprecision(16);
}

}} // namespace precice, io
