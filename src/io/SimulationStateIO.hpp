// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_IO_SIMULATIONSTATEIO_HPP_
#define PRECICE_IO_SIMULATIONSTATEIO_HPP_

#include "logging/Logger.hpp"
#include <string>

namespace precice {
namespace io {

/**
 * @brief Writes global values of the simulation state to a txt file.
 */
class SimulationStateIO
{
public:

  /**
   * @brief Returns standard filename for simulation state file.
   */
  static const std::string& standardFileName();

  /**
   * @brief Returns standard file extension for simulation state file.
   */
  static const std::string& standardFileExtension();

  SimulationStateIO ( const std::string& file );

  void writeState (
    double   time,
    int      timestep,
    long int numberAdvanceCalls );

  void readState (
    double&   time,
    int&      timestep,
    long int& numberAdvanceCalls );

private:

  static logging::Logger _log;

  std::string _file;
};

}} // namespace precice, io

#endif /* PRECICE_IO_SIMULATIONSTATEIO_HPP_ */
