#pragma once

#include "BaseCouplingScheme.hpp"
#include "tarch/logging/Log.h"

namespace precice { namespace cplscheme { namespace tests {
class SerialImplicitCouplingSchemeTest;
} } }


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
    double                      maxTime,
    int                         maxTimesteps,
    double                      timestepLength,
    int                         validDigits,
    const std::string&          firstParticipant,
    const std::string&          secondParticipant,
    const std::string&          localParticipant,
    m2n::M2N::SharedPointer                 m2n,
    constants::TimesteppingMethod dtMethod,
    CouplingMode                cplMode,
    int                         maxIterations = 1
    );

  virtual void initialize(double startTime, int startTimestep);

  virtual void initializeData();

  virtual void advance();

  /// @brief Logging device.
  static tarch::logging::Log _log;

  friend class tests::SerialImplicitCouplingSchemeTest;

};

}}
