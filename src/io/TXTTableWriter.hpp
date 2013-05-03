// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_IO_TXTTABLEWRITER_HPP_
#define PRECICE_IO_TXTTABLEWRITER_HPP_

#include "utils/Globals.hpp"
#include "utils/Dimensions.hpp"
#include <string>
#include <vector>
#include <fstream>

namespace precice {
namespace io {

/**
 * @brief File writer for table-data in text-format.
 *
 * Usage:
 * Create the writer, add data entries in the wanted sequence, and write data
 * values cyclically in the same sequence.
 */
class TXTTableWriter
{
public:

  /**
   * @brief Constants defining possible data types to be written.
   */
  enum DataType {
    INT,
    DOUBLE,
    VECTOR2D,
    VECTOR3D
  };

  /**
   * @brief Constructor, opens file.
   */
  TXTTableWriter ( const std::string& filename );

  /**
   * @brief Destructor, closes file, if not done yet.
   */
  ~TXTTableWriter();

  /**
   * @brief Adds a data entry to the table.
   *
   * Depending on the type and dimension, several text-columns might be added.
   * The writeData() method has to be called in the order in which data entries
   * are added.
   */
  void addData (
    const std::string& name,
    DataType           type );

  /**
   * @brief Writes a integral scalar data value associated to the entry name.
   *
   * The write order is fixed by the order addData() is called.
   */
  void writeData (
    const std::string& name,
    int                value );

  /**
   * @brief Writes a floating-point scalar data value associated to the entry name.
   *
   * The write order is fixed by the order addData() is called.
   */
  void writeData (
    const std::string& name,
    double             value );

  /**
   * @brief Writes a vector data value associated to the entry name.
   *
   * The write order is fixed by the order addData() is called.
   */
  void writeData (
    const std::string&     name,
    const utils::Vector2D& value );

  void writeData (
    const std::string&     name,
    const utils::Vector3D& value );

  /**
   * @brief Closes the file, is automatically called on destruction.
   */
  void close();

private:

  /**
   * @brief Represents one data entry to be written.
   */
  struct Data {

    std::string name;

    DataType type;

    bool operator== ( const Data& data ) const {
      return name == data.name;
    }
  };

  // @brief Logging device.
  static tarch::logging::Log _log;

  std::vector<Data> _data;

  std::vector<Data>::const_iterator _writeIterator;

  std::ofstream _outputStream;

  bool _startedWriting;
};

}} // namespace precice, io

#endif /* PRECICE_IO_TXTTABLEWRITER_HPP_ */
