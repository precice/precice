#ifndef PRECICE_IO_TXTWRITER_HPP_
#define PRECICE_IO_TXTWRITER_HPP_

#include "logging/Logger.hpp"
#include "tarch/la/traits/IsMatrix.h"
#include "tarch/la/traits/IsVector.h"
#include "tarch/la/traits/MatrixTraits.h"
#include "tarch/la/traits/VectorTraits.h"
#include <Eigen/Core>
#include <string>
#include <fstream>

namespace tarch {
  namespace la {
    template<typename Scalar> class DynamicMatrix;
    template<typename Scalar> class DynamicVector;
  }
}

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace io {

/**
 * @brief File writer for matrix in Matlab V7 ASCII format.
 */
class TXTWriter
{
public:

  /**
   * @brief Constructor, opens file and sets format.
   */
  explicit TXTWriter(const std::string& filename);

  /**
   * @brief Destructor, closes file.
   */
  ~TXTWriter();

  ///Writes (appends) the matrix to the file.
  void write(const Eigen::MatrixXd& matrix)
  {
    for (long i = 0; i < matrix.rows(); i++) {
      for (long j = 0; j < matrix.cols(); j++) {
        _file << matrix(i, j) << " ";
      }
    }
    _file << std::endl;
  }

  
private:

  // @brief Logging device.
  static logging::Logger _log;

  // @brief Filestream.
  std::ofstream _file;
};

}} // namespace precice, io

#endif /* PRECICE_IO_TXTWRITER_HPP_ */
