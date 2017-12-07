#pragma once

#include "precice/config/SolverInterfaceConfiguration.hpp"
#include "logging/config/LogConfiguration.hpp"
#include "xml/XMLTag.hpp"
#include "logging/Logger.hpp"

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
   * @brief Returns solver interface configuration.
   */
  const SolverInterfaceConfiguration& getSolverInterfaceConfiguration() const;

private:

  // @brief Logging device.
  static logging::Logger _log;

  // @brief Root tag of preCICE configuration.
  utils::XMLTag _tag;

  LogConfiguration _logConfig;

  SolverInterfaceConfiguration _solverInterfaceConfig;

};

}} // namespace precice, config
