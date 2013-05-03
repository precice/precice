// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef CONFIGURATION_HPP_
#define CONFIGURATION_HPP_

#include "precice/config/SolverInterfaceConfiguration.hpp"
#include "precice/config/LogFilterConfiguration.hpp"
#include "precice/config/LogOutputFormatConfiguration.hpp"
#include "utils/xml/XMLTag.hpp"
#include "tarch/logging/Log.h"

namespace precice {
namespace config {

/**
 * @brief Main class for XML configuration of preCICE.
 */
class Configuration : public utils::XMLTag::Listener
{
public:

  /**
   * @brief Returns the name of the enclosing XML-tag.
   */
  //static const std::string& getTag();

  /**
   * @brief Constructor.
   */
  Configuration();

  /**
   * @brief Returns root xml tag representation.
   */
  utils::XMLTag& getXMLTag();

  /**
   * @brief Reads the information parsed from an xml-file.
   */
  //bool parseSubtag ( utils::XMLTag::XMLReader* xmlReader );

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
   * @brief Returns true, if configuration has validly taken place.
   */
  //bool isValid() const;

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

  utils::XMLTag _tag;

  // @brief Flag to signal validity of configuration.
  //bool _isValid;

  LogFilterConfiguration _logFilterConfig;

  LogOutputFormatConfiguration _logFormatConfig;

  SolverInterfaceConfiguration _solverInterfaceConfig;
};

}} // namespace precice, config

#endif /* CONFIGURATION_HPP_ */
