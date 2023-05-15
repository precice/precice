#pragma once

#include <string>
#include "logging/Logger.hpp"
#include "logging/config/LogConfiguration.hpp"
#include "precice/config/SolverInterfaceConfiguration.hpp"
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

  /**
   * @brief Returns solver interface configuration.
   */
  const SolverInterfaceConfiguration &getSolverInterfaceConfiguration() const;

private:
  logging::Logger _log{"config::Configuration"};

  // @brief Root tag of preCICE configuration.
  xml::XMLTag _tag;

  // The log configuration must be constructed first to prevent log clutter
  LogConfiguration _logConfig;

  // Handle other configuration afterwards
  precice::profiling::ProfilingConfiguration _profilingConfig;

  SolverInterfaceConfiguration _solverInterfaceConfig;
};

} // namespace config
} // namespace precice
