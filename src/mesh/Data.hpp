#pragma once

#include <Eigen/Core>
#include <stddef.h>
#include <string>

#include "logging/Logger.hpp"
#include "precice/types.hpp"

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
   * @brief Do not use this constructor! Only there for compatibility with std::map.
   */
  Data();

  /**
   * @brief Constructor
   */
  Data(
      std::string name,
      DataID      id,
      int         dimension,
      int         spacialDimensions = -1,
      bool        hasGradient       = false);

  /// Returns a reference to the data values.
  Eigen::VectorXd &values();

  /// Returns a const reference to the data values.
  const Eigen::VectorXd &values() const;

  /// Returns a reference to the gradient data values.
  Eigen::MatrixXd &gradientValues();

  /// Returns a const reference to the gradient data values.
  const Eigen::MatrixXd &gradientValues() const;

  /// Returns the name of the data set, as set in the config file.
  const std::string &getName() const;

  /// Returns the ID of the data set (supposed to be unique).
  DataID getID() const;

  /// Sets all values to zero
  void toZero();

  /// Returns if the data contains gradient data
  bool hasGradient() const;

  /// Returns the mesh dimension (i.e., number of rows) of one gradient data value .
  int getSpacialDimensions() const;

  /// Returns the dimension (i.e., number of components) of one data value (i.e number of columns of one gradient data value).
  int getDimensions() const;

private:
  logging::Logger _log{"mesh::Data"};

  Eigen::VectorXd _values;

  Eigen::MatrixXd _gradientValues;

  /// Name of the data set.
  std::string _name;

  /// ID of the data set (supposed to be unique).
  DataID _id;

  /// Spacial Dimension of one element -> number of rows (only 2, 3 allowed for 2D, 3D).
  int _spacialDimensions;

  /// Dimensionality of one data value.
  int _dimensions;

  /// Flag if the gradient data is available
  bool _hasGradient;
};

} // namespace mesh
} // namespace precice
