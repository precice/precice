#include "TXTWriter.hpp"

#include <iomanip>

namespace precice {
namespace io {

TXTWriter:: TXTWriter
(
  const std::string& filename )
:
  _file()
{
  _file.open(filename.c_str());
  if (not _file){
    ERROR("Could not open file \"" << filename << "\" for txt writing!");
  }
  _file.setf ( std::ios::showpoint );
  _file.setf ( std::ios::fixed );
  _file << std::setprecision(16);
}

void TXTWriter::flush()
{
    _file.flush();
}

void TXTWriter::write(const Eigen::MatrixXd& matrix)
{
  for (long i = 0; i < matrix.rows(); i++) {
    for (long j = 0; j < matrix.cols(); j++) {
      _file << matrix(i, j) << ' ';
    }
  }
  _file << '\n';
}

}} // namespace precice, io
