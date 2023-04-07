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

  /// Returns a reference to the data values.
  Eigen::VectorXd &values();

  /// Returns a const reference to the data values.
  const Eigen::VectorXd &values() const;

  /// Returns the name of the data set, as set in the config file.
  const std::string &getName() const;

  /// Returns the ID of the data set (supposed to be unique).
  DataID getID() const;

  /// Sets all values to zero
  void toZero();

  /// Returns the dimension (i.e., number of components) of one data value (i.e number of columns of one gradient data value).
  int getDimensions() const;

private:
  logging::Logger _log{"GlobalData"};

  Eigen::VectorXd _values;

  /// Name of the data set.
  std::string _name;

  /// ID of the data set (supposed to be unique).
  DataID _id;

  /// Dimensionality of one data value.
  int _dimensions;

  /// Spatial Dimension of one element -> number of rows (only 2, 3 allowed for 2D, 3D).
  int _spatialDimensions;
};

} // namespace mesh
} // namespace precice
