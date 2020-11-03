#pragma once

#include <string>
#include <vector>
#include "logging/Logger.hpp"
#include "mesh/Gradient.hpp"
#include "xml/XMLTag.hpp"

namespace precice {
namespace mesh {

/// Performs and provides configuration for Gradient objects from XML files.
class GradientConfiguration : public xml::XMLTag::Listener {
public:
  struct ConfiguredGradient {
    std::string name;
    int         dimensions;

    ConfiguredGradient(
        const std::string &name,
        int                dimensions)
        : name(name), dimensions(dimensions) {}
  };

  GradientConfiguration(xml::XMLTag &parent);

  void setDimensions(int dimensions);

  const std::vector<ConfiguredGradient> &gradients() const;

  ConfiguredGradient getRecentlyConfiguredGradient() const;

  virtual void xmlTagCallback(
      const xml::ConfigurationContext &context,
      xml::XMLTag &                    callingTag);

  virtual void xmlEndTagCallback(
      const xml::ConfigurationContext &context,
      xml::XMLTag &                    callingTag);

  /**
   * @brief Adds gradient manually.
   *
   * @param[in] name Unqiue name of the gradient.
   * @param[in] dataDimensions Dimensionality (1: scalar, 2,3: vector) of data.
   */
  void addGradient(
      const std::string &name,
      int                dataDimensions);

private:
  mutable logging::Logger _log{"mesh::GradientConfiguration"};

  const std::string TAG          = "gradient";
  const std::string ATTR_NAME    = "name";
  const std::string VALUE_VECTOR = "vector";
  const std::string VALUE_SCALAR = "scalar";

  /// Dimension of space.
  int _dimensions = 0;

  std::vector<ConfiguredGradient> _gradients;

  int _indexLastConfigured = -1;

  int getDataDimensions(const std::string &typeName) const;
};

} // namespace mesh
} // namespace precice
