#pragma once

#include <Eigen/Core>
#include <fstream>
#include <string>
#include "logging/Logger.hpp"

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace io {

/**
 * @brief File reader for matrix/vector in Matlab V7 ASCII format.
 */
class TXTReader {
public:
  /// Constructor, opens file and sets format.
  explicit TXTReader(const std::string &filename);

  /// Reads the Eigen::Matrix from the file.
  template <typename Scalar, int Rows, int Cols>
  void read(Eigen::Matrix<Scalar, Rows, Cols> &matrix)
  {
    for (long i = 0; i < matrix.rows(); i++) {
      for (long j = 0; j < matrix.cols(); j++) {
        double scalar;
        _file >> scalar;
        matrix(i, j) = scalar;
      }
    }
  }

private:
  logging::Logger _log{"io::TXTReader"};

  /// @brief Filestream.
  std::ifstream _file;
};

} // namespace io
} // namespace precice
