#pragma once

#include <string>
#include <vector>
#include "logging/Logger.hpp"
#include "mesh/Data.hpp"
#include "time/Time.hpp"
#include "mesh/SharedPointer.hpp"
#include "utils/ManageUniqueIDs.hpp"
#include "xml/XMLTag.hpp"

namespace precice {
namespace mesh {

/// Performs and provides configuration for Data objects from XML files.
class DataConfiguration : public xml::XMLTag::Listener {
public:
  struct ConfiguredData {
    std::string name;
    int         dimensions;
    int         waveformDegree;
    bool        isGlobal; // false = mesh-associated data, true = global data

    ConfiguredData(
        const std::string &name,
        int                dimensions,
        int                waveformDegree,
        bool               isGlobal)
        : name(name), dimensions(dimensions), waveformDegree(waveformDegree), isGlobal(isGlobal) {}
  };

  DataConfiguration(xml::XMLTag &parent);

  void setDimensions(int dimensions);

  const std::vector<ConfiguredData> &data() const;

  /// Returns the Global Data with the matching name
  const PtrData &globalData(const std::string &dataName) const;

  /// Returns whether Data Configuration has Global Data with the dataName
  bool hasGlobalDataName(const std::string &dataName) const;

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
  void addData(const std::string &name,
               int                dataDimensions,
               int                waveformDegree = time::Time::DEFAULT_WAVEFORM_DEGREE,
               bool               isGlobal);

  void setExperimental(bool experimental);

private:
  mutable logging::Logger _log{"mesh::DataConfiguration"};

  const std::string TAG_MESH_DATA   = "data";
  const std::string TAG_GLOBAL_DATA = "global-data";
  const std::string ATTR_NAME    = "name";
  const std::string ATTR_DEGREE  = "waveform-degree";
  const std::string VALUE_VECTOR = "vector";
  const std::string VALUE_SCALAR = "scalar";

  /// Dimension of space.
  int _dimensions = 0;

  std::vector<ConfiguredData> _data;

  std::vector<PtrData> _globalData;

  int _indexLastConfigured = -1;

  utils::ManageUniqueIDs _dataIDManager;

  int getDataDimensions(const std::string &typeName) const;

  /// Creates a GlobalData object and appends it to _globalData vector
  void createGlobalData(const std::string &name,
                        int                dimension,
                        DataID             id);

  bool _experimental = false;
};

} // namespace mesh
} // namespace precice
