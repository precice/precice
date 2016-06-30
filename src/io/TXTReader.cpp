#include "TXTReader.hpp"
#include "utils/Globals.hpp"
#include "tarch/la/DynamicMatrix.h"
#include "tarch/la/DynamicVector.h"

namespace precice {
namespace io {

tarch::logging::Log TXTReader:: _log ( "precice::io::TXTReader" );

TXTReader:: TXTReader
(
  const std::string& filename )
:
  _file()
{
  _file.open(filename.c_str());
  if (not _file){
    preciceError("TXTReader()", "Could not open file \"" << filename
                 << "\" for txt reading!");
  }
  _file.setf(std::ios::showpoint);
  _file.setf(std::ios::fixed);
  //_file << std::setprecision(16);
}

TXTReader:: ~TXTReader()
{
  if (_file){
    _file.close();
  }
}




}} // namespace precice, io
