#pragma once

#include <string>
#include <vector>
#include "logging/Logger.hpp"
#include "mesh/Data.hpp"
#include "time/Time.hpp"
#include "xml/XMLTag.hpp"

namespace precice {
namespace mesh {

/// Performs and provides configuration for Data objects from XML files.
class DataConfiguration : public xml::XMLTag::Listener {
public:
  struct ConfiguredData {
    std::string name;
    int         dimensions;
    int         interpolationOrder;

    ConfiguredData(
        const std::string &name,
        int                dimensions,
        int                interpolationOrder = time::Time::UNDEFINED_INTERPOLATION_ORDER)
        : name(name), dimensions(dimensions), interpolationOrder(interpolationOrder) {}
  };

  DataConfiguration(xml::XMLTag &parent);

  void setDimensions(int dimensions);

  const std::vector<ConfiguredData> &data() const;

  ConfiguredData getRecentlyConfiguredData() const;

  virtual void xmlTagCallback(
      const xml::ConfigurationContext &context,
      xml::XMLTag &                    callingTag);

  virtual void xmlEndTagCallback(
      const xml::ConfigurationContext &context,
      xml::XMLTag &                    callingTag);

  /**
   * @brief Adds data manually.
   *
   * @param[in] name Unqiue name of the data.
   * @param[in] dataDimensions Dimensionality (1: scalar, 2,3: vector) of data.
   */
  void addData(const std::string &name,
               int                dataDimensions,
               int                interpolationOrder = time::Time::UNDEFINED_INTERPOLATION_ORDER);

private:
  mutable logging::Logger _log{"mesh::DataConfiguration"};

  const std::string TAG          = "data";
  const std::string ATTR_NAME    = "name";
  const std::string TAG_WAVEFORM = "waveform";
  const std::string ATTR_ORDER   = "order";
  const std::string VALUE_VECTOR = "vector";
  const std::string VALUE_SCALAR = "scalar";

  struct Config {
    std::string name;
    int         dataDimensions;
    int         waveformOrder = time::Time::UNDEFINED_INTERPOLATION_ORDER;
  } _config;

  /// Dimension of space.
  int _dimensions = 0;

  std::vector<ConfiguredData> _data;

  int _indexLastConfigured = -1;

  int getDataDimensions(const std::string &typeName) const;
};

} // namespace mesh
} // namespace precice
