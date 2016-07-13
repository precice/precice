#pragma once

#include "precice/config/SolverInterfaceConfiguration.hpp"
#include "precice/config/LogFilterConfiguration.hpp"
#include "precice/config/LogOutputFormatConfiguration.hpp"
#include "utils/xml/XMLTag.hpp"
#include "tarch/logging/Log.h"

namespace precice {
namespace config {

/**
 * @brief Main class for preCICE XML configuration tree.
 *
 * The configuration process is triggered by fetching the root tag with method
 * getXMLTag() and calling its parse() method.
 */
class Configuration : public utils::XMLTag::Listener
{
public:

  /**
   * @brief Constructor.
   */
  Configuration();

  /**
   * @brief Destructor, empty.
   */
  virtual ~Configuration() {}

  /**
   * @brief Returns root xml tag to start the automatic configuration process.
   */
  utils::XMLTag& getXMLTag();

  /**
   * @brief Callback function required for use of automatic configuration.
   *
   * @return True, if successful.
   */
  virtual void xmlTagCallback ( utils::XMLTag& tag );

  /**
   * @brief Callback function required for use of automatic configuration.
   *
   * @return True, if successful.
   */
  virtual void xmlEndTagCallback ( utils::XMLTag& tag );

  /**
   * @brief Returns log filter configuration.
   */
  const LogFilterConfiguration& getLogFilterConfiguration() const;

  /**
   * @brief Returns log output format configuration.
   */
  const LogOutputFormatConfiguration& getLogFormatConfiguration() const;

  /**
   * @brief Returns solver interface configuration.
   */
  const SolverInterfaceConfiguration& getSolverInterfaceConfiguration() const;

private:

  // @brief Logging device.
  static tarch::logging::Log _log;

  // @brief Root tag of preCICE configuration.
  utils::XMLTag _tag;

  LogFilterConfiguration _logFilterConfig;

  LogOutputFormatConfiguration _logFormatConfig;

  SolverInterfaceConfiguration _solverInterfaceConfig;
};

}} // namespace precice, config
