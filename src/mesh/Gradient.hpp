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
 * @brief Describes a set of gradient values belonging to the vertices of a mesh.
 */
class Gradient {
public:

  /**
   * @brief Returns the number of created (and still existing) Gradient objects.
   *
   * Used to give Gradient objects unique IDs.
   */
  static size_t getGradientCount();

  /**
   * @brief Sets the gradient counter to zero.
   *
   * Used in test cases where multiple scenarios with gradients are run.
   */
  static void resetGradientCount();

  /**
   * @brief Do not use this constructor! Only there for compatibility with std::map.
   */
  Gradient();

  /**
   * @brief Constructor.
   */
  Gradient(
      const std::string &name,
      int                id,
      int                dimension);

  /// Destructor, decrements gradient count.
  ~Gradient();

  /// Returns a reference to the gradient values.
  Eigen::MatrixXd &values();

  /// Returns a const reference to the gradient values.
  const Eigen::MatrixXd &values() const;

  /// Returns the name of the gradient data set, as set in the config file.
  const std::string &getName() const;

  /// Returns the ID of the gradient data set (supposed to be unique).
  int getID() const;

  /// Sets all values to zero
  void toZero();

  /// Returns the dimension (i.e., number of components) of one gradient value.
  int getDimensions() const;

private:
  logging::Logger _log{"mesh::Gradient"};

  /// Counter for existing Gradient objects.
  static size_t _gradientCount;

  /// The gradients matrices. Column-major ordering is default.
  Eigen::MatrixXd _values;

  /// Name of the gradient data set.
  std::string _name;

  /// ID of the gradient data set (supposed to be unique).
  int _id;

  /// Dimensionality of one gradient value.
  int _dimensions;

};

} // namespace mesh
} // namespace precice
