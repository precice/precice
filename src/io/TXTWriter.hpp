#pragma once

#include "logging/Logger.hpp"
#include <Eigen/Core>
#include <string>
#include <fstream>

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
  void write(const Eigen::MatrixXd& matrix);

  ///Flush the buffer to file
  void flush();

  
private:
  logging::Logger _log{"io::TXTWriter"};

  // @brief Filestream.
  std::ofstream _file;
};

}} // namespace precice, io

