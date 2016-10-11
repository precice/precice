#ifndef PRECICE_IO_TXTREADER_HPP_
#define PRECICE_IO_TXTREADER_HPP_

#include "logging/Logger.hpp"
#include "tarch/la/traits/IsMatrix.h"
#include "tarch/la/traits/IsVector.h"
#include "tarch/la/traits/MatrixTraits.h"
#include "tarch/la/traits/VectorTraits.h"
#include <Eigen/Core>
#include <string>
#include <fstream>

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace io {

/**
 * @brief File reader for matrix/vector in Matlab V7 ASCII format.
 */
class TXTReader
{
public:

  /**
   * @brief Constructor, opens file and sets format.
   */
  explicit TXTReader(const std::string& filename);

  /**
   * @brief Destructor, closes file.
   */
  ~TXTReader();

  /// Reads the matrix from the file.
  template<typename MATRIX>
  typename std::enable_if<tarch::la::IsMatrix<MATRIX>::value
  >::type read(MATRIX& matrix)
  {
    typedef tarch::la::MatrixTraits<MATRIX> T;
    for (int i=0; i < T::rows(matrix); i++){
      for (int j=0; j < T::cols(matrix); j++){
        _file >> T::elem(i,j,matrix);
      }
    }
  }

  /// Reads the vector from the file.
  template<typename VECTOR>
  typename std::enable_if<tarch::la::IsVector<VECTOR>::value
  >::type read(VECTOR& vector)
  {
    typedef tarch::la::VectorTraits<VECTOR> T;
    for (int i=0; i < T::size(vector); i++){
      _file >> T::elem(i,vector);
    }
  }

  /// Reads the Eigen::Matrix from the file.
  template<typename Scalar, int Rows, int Cols>
  void read(Eigen::Matrix<Scalar, Rows, Cols>& matrix)
  {
    for (long i = 0; i < matrix.rows(); i++) {
      for (long j = 0; j < matrix.cols(); j++) {
        double scalar;
        _file >> scalar;
        matrix(i,j) = scalar;
      }
    }
  }


private:

  /// @brief Logging device.
  static logging::Logger _log;

  /// @brief Filestream.
  std::ifstream _file;
};

}} // namespace precice, io

#endif /* PRECICE_IO_TXTREADER_HPP_ */
