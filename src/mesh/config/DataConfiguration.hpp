#pragma once

#include <string>
#include <vector>
#include "logging/Logger.hpp"
#include "mesh/Data.hpp"
#include "time/Time.hpp"
#include "utils/ManageUniqueIDs.hpp"
#include "xml/XMLTag.hpp"

namespace precice::mesh {

/// Performs and provides configuration for Data objects from XML files.
class DataConfiguration : public xml::XMLTag::Listener {
public:
  struct ConfiguredData {
    std::string                        name;
    Data::typeName                     typeName;
    int                                waveformDegree;
    std::vector<std::optional<double>> lowerBound;
    std::vector<std::optional<double>> upperBound;

    ConfiguredData(
        const std::string                 &name,
        const Data::typeName               typeName,
        int                                waveformDegree,
        std::vector<std::optional<double>> lowerBound,
        std::vector<std::optional<double>> upperBound)
        : name(name), typeName(typeName), waveformDegree(waveformDegree), lowerBound(lowerBound), upperBound(upperBound) {}
  };

  DataConfiguration(xml::XMLTag &parent);

  const std::vector<ConfiguredData> &data() const;

  ConfiguredData getRecentlyConfiguredData() const;

  void xmlTagCallback(
      const xml::ConfigurationContext &context,
      xml::XMLTag                     &callingTag) override;

  void xmlEndTagCallback(
      const xml::ConfigurationContext &context,
      xml::XMLTag                     &callingTag) override;

  /**
   * @brief Adds data manually.
   *
   * @param[in] name Unique name of the data.
   * @param[in] dataDimensions Dimensionality (1: scalar, 2,3: vector) of data.
   * @param[in] waveformDegree Degree of waveform associated with this data.
   * @param[in] lowerBound Lower bound of the data for violation check.
   * @param[in] upperBound Upper bound of the data for violation check.
   */
  void addData(const std::string                 &name,
               const Data::typeName               typeName,
               int                                waveformDegree = time::Time::DEFAULT_WAVEFORM_DEGREE,
               std::vector<std::optional<double>> lowerBound     = std::vector<std::optional<double>>(3),
               std::vector<std::optional<double>> upperBound     = std::vector<std::optional<double>>(3));

private:
  mutable logging::Logger _log{"mesh::DataConfiguration"};

  const std::string TAG                = "data";
  const std::string ATTR_NAME          = "name";
  const std::string ATTR_DEGREE        = "waveform-degree";
  const std::string ATTR_LOWER_BOUND   = "lower-bound";
  const std::string ATTR_UPPER_BOUND   = "upper-bound";
  const std::string ATTR_LOWER_BOUND_X = "lower-bound-x";
  const std::string ATTR_LOWER_BOUND_Y = "lower-bound-y";
  const std::string ATTR_LOWER_BOUND_Z = "lower-bound-z";
  const std::string ATTR_UPPER_BOUND_X = "upper-bound-x";
  const std::string ATTR_UPPER_BOUND_Y = "upper-bound-y";
  const std::string ATTR_UPPER_BOUND_Z = "upper-bound-z";
  const std::string VALUE_VECTOR       = "vector";
  const std::string VALUE_SCALAR       = "scalar";

  std::vector<ConfiguredData> _data;

  int _indexLastConfigured = -1;
};

} // namespace precice::mesh
