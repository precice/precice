#pragma once

#include "mesh/Data.hpp"
#include "xml/XMLTag.hpp"
#include "logging/Logger.hpp"
#include <vector>
#include <string>

namespace precice {
namespace mesh {

/// Performs and provides configuration for Data objects from XML files.
class DataConfiguration : public xml::XMLTag::Listener
{
public:

  struct ConfiguredData
  {
    std::string name;
    int dimensions;

    ConfiguredData (
      const std::string& name,
      int                dimensions )
    : name(name), dimensions(dimensions) {}
  };

  /**
   * @brief Returns the name of the main XML-tag of this configuration.
   */
  //static const std::string& getTag();

  DataConfiguration ( xml::XMLTag& parent );

  void setDimensions ( int dimensions );

  /**
   * @brief Returns true, if configuration was successfull.
   */
  //bool isValid() const;

  const std::vector<ConfiguredData>& data() const;

  ConfiguredData getRecentlyConfiguredData() const;

  virtual void xmlTagCallback ( xml::XMLTag& callingTag );

  virtual void xmlEndTagCallback ( xml::XMLTag& callingTag );

  /**
   * @brief Adds data manually.
   *
   * @param[in] name Unqiue name of the data.
   * @param[in] dataDimensions Dimensionality (1: scalar, 2,3: vector) of data.
   */
  void addData (
    const std::string& name,
    int                dataDimensions );

  //int getDimensions() const;

private:
  
  mutable logging::Logger _log{"mesh::DataConfiguration"};

  const std::string TAG;
  const std::string ATTR_NAME;
  const std::string VALUE_VECTOR;
  const std::string VALUE_SCALAR;

  /// Dimension of space.
  int _dimensions;

  std::vector<ConfiguredData> _data;

  int _indexLastConfigured;

  int getDataDimensions(const std::string& typeName) const;
};

}} // namespace precice, mesh
