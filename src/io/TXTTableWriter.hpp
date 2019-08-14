#pragma once

#include <Eigen/Core>
#include <fstream>
#include <string>
#include <vector>
#include "logging/Logger.hpp"

namespace precice {
namespace io {

/**
 * @brief File writer for table-data in text-format.
 *
 * Usage:
 * Create the writer, add data entries in the wanted sequence, and write data
 * values cyclically in the same sequence.
 */
class TXTTableWriter {
public:
  /// Constants defining possible data types to be written.
  enum DataType {
    INT,
    DOUBLE,
    VECTOR2D,
    VECTOR3D
  };

  /// Constructor, opens file.
  explicit TXTTableWriter(const std::string &filename);

  /**
   * @brief Adds a data entry to the table.
   *
   * Depending on the type and dimension, several text-columns might be added.
   * The writeData() method has to be called in the order in which data entries
   * are added.
   */
  void addData(
      const std::string &name,
      DataType           type);

  /**
   * @brief Writes a integral scalar data value associated to the entry name.
   *
   * The write order is fixed by the order addData() is called.
   */
  void writeData(
      const std::string &name,
      int                value);

  /**
   * @brief Writes a floating-point scalar data value associated to the entry name.
   *
   * The write order is fixed by the order addData() is called.
   */
  void writeData(
      const std::string &name,
      double             value);

  /**
   * @brief Writes a vector data value associated to the entry name.
   *
   * The write order is fixed by the order addData() is called.
   */
  void writeData(
      const std::string &    name,
      const Eigen::Vector2d &value);

  void writeData(
      const std::string &    name,
      const Eigen::Vector3d &value);

  /// Closes the file, is automatically called on destruction.
  void close();

  /// Resets the table information.
  void reset();

private:
  /// Represents one data entry to be written.
  struct Data {

    std::string name;

    DataType type;

    bool operator==(const Data &data) const
    {
      return name == data.name;
    }
  };

  logging::Logger _log{"io::TXTTableWriter"};

  std::vector<Data> _data;

  std::vector<Data>::const_iterator _writeIterator;

  std::ofstream _outputStream;
};

} // namespace io
} // namespace precice
