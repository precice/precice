// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
extern "C" {
#include "Constants.h"
}

#include "precice/Constants.hpp"

const char* precicec_nameConfiguration()
{
  return precice::constants::nameConfiguration().c_str();
}

const char* precicec_actionWriteIterationCheckpoint()
{
  return precice::constants::actionWriteIterationCheckpoint().c_str();
}

const char* precicec_actionReadIterationCheckpoint()
{
  return precice::constants::actionReadIterationCheckpoint().c_str();
}

const char* precicec_actionWriteSimulationCheckpoint()
{
  return precice::constants::actionWriteSimulationCheckpoint().c_str();
}

const char* precicec_actionReadSimulationCheckpoint()
{
  return precice::constants::actionReadSimulationCheckpoint().c_str();
}

//int precice_dimensions = precice::constants::dimensions();
