#pragma once

#include <Eigen/Core>
#include <stddef.h>
#include <string>

#include "logging/Logger.hpp"
#include "precice/types.hpp"

namespace precice {
namespace mesh {
class GlobalData;
}
} // namespace precice
// TODO: GlobalData should not be inside mesh namespace

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace mesh {
/**
 * @brief Describes a data value not belonging to a mesh.
 */
class GlobalData {
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
  GlobalData();

  /**
   * @brief Constructor
   */
  GlobalData(
      std::string name,
      DataID      id,
      int         dimension,
      int         spatialDimensions = -1);

  /// Returns a const reference to the data values.
  const Eigen::VectorXd &values() const;

//   no gradients for meshless data
//   Eigen::MatrixXd &gradientValues();

//   no gradients for meshless data
//   const Eigen::MatrixXd &gradientValues() const;

  /// Returns the name of the data set, as set in the config file.
  const std::string &getName() const;

  /// Returns the ID of the data set (supposed to be unique).
  DataID getID() const;

  /// Sets all values to zero
  void toZero();

//   no gradients for meshless data
//   bool hasGradient() const;

//   no gradients for meshless data
//   void requireDataGradient();

  /// Returns the mesh dimension (i.e., number of rows) of one gradient data value .
  int getSpatialDimensions() const;

  /// Returns the dimension (i.e., number of components) of one data value (i.e number of columns of one gradient data value).
  int getDimensions() const;


  private:
  logging::Logger _log{"GlobalData"};

  Eigen::VectorXd _values;

//   no gradients for meshless data
//   Eigen::MatrixXd _gradientValues;

  /// Name of the data set.
  std::string _name;

  /// ID of the data set (supposed to be unique).
  DataID _id;

  /// Dimensionality of one data value.
  int _dimensions;

  /// Spatial Dimension of one element -> number of rows (only 2, 3 allowed for 2D, 3D).
  int _spatialDimensions;

//   no gradients for meshless data
//   bool _hasGradient = false;

};

} // namespace mesh
} // namespace precice