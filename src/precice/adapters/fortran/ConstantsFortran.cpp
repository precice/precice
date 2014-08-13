/* Copyright (C) 2011 Technische Universitaet Muenchen
 * This file is part of the preCICE project. For conditions of distribution and
 * use, please see the license notice at http://www5.in.tum.de/wiki/index.php/precice_c_License */

extern"C" {
#include "ConstantsFortran.hpp"
}
#include "precice/Constants.hpp"
#include "tarch/Assertions.h"

void precicef_name_config_
(
  char*  nameConfig,
  int lengthNameConfig )
{
  const std::string& name = precice::constants::nameConfiguration();
  assertion2(name.size() < lengthNameConfig, name.size(), lengthNameConfig);
  for (size_t i=0; i < name.size(); i++){
    nameConfig[i] = name[i];
  }
}

void precicef_action_write_iter_checkp_
(
  char*  nameAction,
  int lengthNameAction )
{
  const std::string& name = precice::constants::actionWriteIterationCheckpoint();
  assertion2(name.size() < lengthNameAction, name.size(), lengthNameAction);
  for (size_t i=0; i < name.size(); i++){
    nameAction[i] = name[i];
  }
}

void precicef_action_write_initial_data_(
  char*  nameAction,
  int lengthNameAction )
{
  const std::string& name = precice::constants::actionWriteInitialData();
  assertion2(name.size() < lengthNameAction, name.size(), lengthNameAction);
  for (size_t i=0; i < name.size(); i++){
    nameAction[i] = name[i];
  }
}

void precicef_action_read_iter_checkp_
(
  char*  nameAction,
  int lengthNameAction )
{
  const std::string& name = precice::constants::actionReadIterationCheckpoint();
  assertion2(name.size() < lengthNameAction, name.size(), lengthNameAction);
  for (size_t i=0; i < name.size(); i++){
    nameAction[i] = name[i];
  }
}


void precicef_action_write_sim_checkp_
(
  char*  nameAction,
  int lengthNameAction )
{
  const std::string& name = precice::constants::actionWriteSimulationCheckpoint();
  assertion2(name.size() < lengthNameAction, name.size(), lengthNameAction);
  for (size_t i=0; i < name.size(); i++){
    nameAction[i] = name[i];
  }
}

void precicef_action_read_sim_checkp_
(
  char*  nameAction,
  int lengthNameAction )
{
  const std::string& name = precice::constants::actionReadSimulationCheckpoint();
  assertion2(name.size() < lengthNameAction, name.size(), lengthNameAction);
  for (size_t i=0; i < name.size(); i++){
    nameAction[i] = name[i];
  }
}
