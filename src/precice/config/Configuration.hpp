#pragma once

#include <string>
#include "SharedPointer.hpp"
#include "boost/smart_ptr.hpp"
#include "cplscheme/SharedPointer.hpp"
#include "logging/Logger.hpp"
#include "logging/config/LogConfiguration.hpp"
#include "m2n/M2N.hpp"
#include "m2n/config/M2NConfiguration.hpp"
#include "mapping/SharedPointer.hpp"
#include "mesh/SharedPointer.hpp"
#include "profiling/config/ProfilingConfiguration.hpp"
#include "xml/XMLTag.hpp"

namespace precice {
namespace config {

/**
 * @brief Main class for preCICE XML configuration tree.
 *
 * The configuration process is triggered by fetching the root tag with method
 * getXMLTag() and calling its parse() method.
 */
class Configuration : public xml::XMLTag::Listener {
public:
  Configuration();

  /**
   * @brief Destructor, empty.
   */
  virtual ~Configuration() {}

  /**
   * @brief Returns root xml tag to start the automatic configuration process.
   */
  xml::XMLTag &getXMLTag();

  /**
   * @brief Callback function required for use of automatic configuration.
   *
   * @return True, if successful.
   */
  virtual void xmlTagCallback(const xml::ConfigurationContext &context, xml::XMLTag &tag);

  /**
   * @brief Callback function required for use of automatic configuration.
   *
   * @return True, if successful.
   */
  virtual void xmlEndTagCallback(const xml::ConfigurationContext &context, xml::XMLTag &tag);

  /// @brief Returns whether experimental features are allowed or not
  bool allowsExperimental() const
  {
    return _experimental;
  }

  /// @brief Returns whether participants wait for each other in finalize
  bool waitInFinalize() const
  {
    return _waitInFinalize;
  }

  const mesh::PtrDataConfiguration getDataConfiguration() const
  {
    return _dataConfiguration;
  }

  const mesh::PtrMeshConfiguration getMeshConfiguration() const
  {
    return _meshConfiguration;
  }

  const m2n::M2NConfiguration::SharedPointer getM2NConfiguration() const
  {
    return _m2nConfiguration;
  }

  const PtrParticipantConfiguration &getParticipantConfiguration() const;

  const cplscheme::PtrCouplingSchemeConfiguration getCouplingSchemeConfiguration() const
  {
    return _couplingSchemeConfiguration;
  }

  /**
   * @brief For manual configuration in test cases.
   */
  void setDataConfiguration(mesh::PtrDataConfiguration config)
  {
    _dataConfiguration = config;
  }

  /**
   * @brief For manual configuration in test cases.
   */
  void setMeshConfiguration(mesh::PtrMeshConfiguration config)
  {
    _meshConfiguration = config;
  }

  /**
    * @brief For manual configuration in test cases.
    */
  void setParticipantConfiguration(PtrParticipantConfiguration config)
  {
    _participantConfiguration = config;
  }

private:
  logging::Logger _log{"config::Configuration"};

  /// Allow the use of experimental features
  bool _experimental = false;

  /// Synchronize participants in finalize
  bool _waitInFinalize;

  // @brief Root tag of preCICE configuration.
  xml::XMLTag _tag;

  // The log configuration must be constructed first to prevent log clutter
  LogConfiguration _logConfig;

  // Handle other configuration afterwards
  precice::profiling::ProfilingConfiguration _profilingConfig;

  mesh::PtrDataConfiguration _dataConfiguration;

  mesh::PtrMeshConfiguration _meshConfiguration;

  m2n::M2NConfiguration::SharedPointer _m2nConfiguration;

  PtrParticipantConfiguration _participantConfiguration;

  cplscheme::PtrCouplingSchemeConfiguration _couplingSchemeConfiguration;
};

} // namespace config
} // namespace precice
