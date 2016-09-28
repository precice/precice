#ifndef PRECICE_MESH_DATACONFIGURATION_HPP_
#define PRECICE_MESH_DATACONFIGURATION_HPP_

#include "mesh/SharedPointer.hpp"
#include "mesh/Data.hpp"
#include "utils/xml/XMLTag.hpp"
#include "logging/Logger.hpp"
#include <vector>
#include <string>

namespace precice {
namespace mesh {

/**
 * @brief Performs and provides configuration for Data objects from XML files.
 */
class DataConfiguration : public utils::XMLTag::Listener
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

  /**
   * @brief Constructor.
   */
  DataConfiguration ( utils::XMLTag& parent );

  void setDimensions ( int dimensions );

  /**
   * @brief Returns true, if configuration was successfull.
   */
  //bool isValid() const;

  const std::vector<ConfiguredData>& data() const;

  ConfiguredData getRecentlyConfiguredData() const;

  virtual void xmlTagCallback ( utils::XMLTag& callingTag );

  virtual void xmlEndTagCallback ( utils::XMLTag& callingTag );

  /**
   * @brief Adds data manually.
   *
   * @param name [IN] Unqiue name of the data.
   * @param dataDimensions [IN] Dimensionality (1: scalar, 2,3: vector) of data.
   */
  void addData (
    const std::string& name,
    int                dataDimensions );

  //int getDimensions() const;

private:

  static logging::Logger _log;

  const std::string TAG;
  const std::string ATTR_NAME;
  //const std::string ATTR_TYPE;
  const std::string VALUE_VECTOR;
  const std::string VALUE_SCALAR;

//  utils::XMLTag _tag;

  // @brief Dimension of space.
  int _dimensions;

  //bool _isValid;

  std::vector<ConfiguredData> _data;

  int _indexLastConfigured;

  int getDataDimensions(const std::string& typeName) const;
};

}} // namespace precice, mesh

#endif /* PRECICE_MESH_DATACONFIGURATION_HPP_ */
