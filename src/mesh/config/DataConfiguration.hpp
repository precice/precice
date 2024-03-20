#pragma once

#include <string>
#include <vector>
#include "logging/Logger.hpp"
#include "mesh/Data.hpp"
#include "time/Time.hpp"
#include "utils/ManageUniqueIDs.hpp"
#include "xml/XMLTag.hpp"

namespace precice {
namespace mesh {

/// Performs and provides configuration for Data objects from XML files.
class DataConfiguration : public xml::XMLTag::Listener {
public:
  struct ConfiguredData {
    std::string    name;
    Data::typeName typeName;
    int            waveformDegree;

    ConfiguredData(
        const std::string &  name,
        const Data::typeName typeName,
        int                  waveformDegree)
        : name(name), typeName(typeName), waveformDegree(waveformDegree) {}
  };

  DataConfiguration(xml::XMLTag &parent);

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
   * @param[in] name Unique name of the data.
   * @param[in] dataDimensions Dimensionality (1: scalar, 2,3: vector) of data.
   * @param[in] waveformDegree Degree of waveform associated with this data.
   */
  void addData(const std::string &  name,
               const Data::typeName typeName,
               int                  waveformDegree = time::Time::DEFAULT_WAVEFORM_DEGREE);

private:
  mutable logging::Logger _log{"mesh::DataConfiguration"};

  const std::string TAG          = "data";
  const std::string ATTR_NAME    = "name";
  const std::string ATTR_DEGREE  = "waveform-degree";
  const std::string VALUE_VECTOR = "vector";
  const std::string VALUE_SCALAR = "scalar";

  std::vector<ConfiguredData> _data;

  int _indexLastConfigured = -1;
};

} // namespace mesh
} // namespace precice
