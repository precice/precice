// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef CONFIGURATION_HPP_
#define CONFIGURATION_HPP_

#include "precice/config/SolverInterfaceConfiguration.hpp"
#include "precice/config/LogFilterConfiguration.hpp"
#include "precice/config/LogOutputFormatConfiguration.hpp"
#include "utils/xml/XMLTag.hpp"
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
  static logging::Logger _log;

  // @brief Root tag of preCICE configuration.
  utils::XMLTag _tag;

  LogFilterConfiguration _logFilterConfig;

  LogOutputFormatConfiguration _logFormatConfig;

  SolverInterfaceConfiguration _solverInterfaceConfig;
};

}} // namespace precice, config

#endif /* CONFIGURATION_HPP_ */
