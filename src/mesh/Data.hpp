#pragma once

#include <Eigen/Core>
#include <stddef.h>
#include <string>
#include "logging/Logger.hpp"

namespace precice {
namespace mesh {
class Mesh;
}
} // namespace precice

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace mesh {

/**
 * @brief Describes a set of data values belonging to the vertices of a mesh.
 */
class Data {
public:
  // @brief Possible types of data values.
  //  enum DataType {
  //    TYPE_UNDEFINED,
  //    TYPE_DOUBLE,
  //    TYPE_VECTOR
  //  };

  // @brief Name of an undefined data type.
  //static const std::string TYPE_NAME_UNDEFINED;
  // @brief Name of a double data type.
  //static const std::string TYPE_NAME_DOUBLE;
  // @brief Name of a vector data type.
  //static const std::string TYPE_NAME_VECTOR;

  /**
   * @brief Returns the number of created (and still existing) Data objects.
   *
   * Used to give Data objects unique IDs.
   */
  static size_t getDataCount();

  /**
   * @brief Sets the data counter to zero.
   *
   * Used in test cases where multiple scenarios with data are run.
   */
  static void resetDataCount();

  /**
   * @brief Do not use this constructor! Only there for compatibility with std::map.
   */
  Data();

  /**
   * @brief Constructor.
   */
  Data(
      const std::string &name,
      int                id,
      int                dimension);

  /// Destructor, decrements data count.
  ~Data();

  /// Returns a reference to the data values.
  Eigen::VectorXd &values();

  /// Returns a const reference to the data values.
  const Eigen::VectorXd &values() const;

  /// Returns the name of the data set, as set in the config file.
  const std::string &getName() const;

  /// Returns the ID of the data set (supposed to be unique).
  int getID() const;

  /**
   * @brief Returns the type constant of the data set.
   */
  //  DataType getType () const;

  /// Sets all values to zero
  void toZero();

  /// Returns the dimension (i.e., number of components) of one data value.
  int getDimensions() const;

private:
  logging::Logger _log{"mesh::Data"};

  /// Counter for existing Data objects.
  static size_t _dataCount;

  Eigen::VectorXd _values;

  /// Name of the data set.
  std::string _name;

  /// ID of the data set (supposed to be unique).
  int _id;

  //  // @brief Type of data (scalar or vector).
  //  DataType _type;

  /// Dimensionality of one data value.
  int _dimensions;
};

} // namespace mesh
} // namespace precice
