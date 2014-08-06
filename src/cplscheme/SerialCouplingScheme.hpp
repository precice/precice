#pragma once

#include "BaseCouplingScheme.hpp"
#include "tarch/logging/Log.h"
#include "com/Communication.hpp"
#include "com/SharedPointer.hpp"


namespace precice {
namespace cplscheme {

class SerialCouplingScheme : public BaseCouplingScheme
{
public:
  
/**
 * @brief Constructor.
 *
 * @param maxTime [IN] Simulation time limit, or UNDEFINED_TIME.
 * @param maxTimesteps [IN] Simulation timestep limit, or UNDEFINED_TIMESTEPS.
 * @param timestepLength [IN] Simulation timestep length.
 * @param firstParticipant [IN] Name of participant starting simulation.
 * @param secondParticipant [IN] Name of second participant in coupling.
 * @param localParticipant [IN] Name of participant using this coupling scheme.
 * @param communication [IN] Communication object for com. between participants.
 * @param monitorIterations [IN] If true, a txt file monitoring iterations is written.
 */
SerialCouplingScheme (
  double                maxTime,
  int                   maxTimesteps,
  double                timestepLength,
  int                   validDigits,
  const std::string&    firstParticipant,
  const std::string&    secondParticipant,
  const std::string&    localParticipant,
  com::PtrCommunication communication,
  int                   maxIterations,
  constants::TimesteppingMethod dtMethod);

/// @brief Logging device.
static tarch::logging::Log _log;

};

}}
