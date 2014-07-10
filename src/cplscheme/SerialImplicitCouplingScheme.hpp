// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_CPLSCHEME_SERIALIMPLICITCOUPLINGSCHEME_HPP_
#define PRECICE_CPLSCHEME_SERIALIMPLICITCOUPLINGSCHEME_HPP_

#include "ImplicitCouplingScheme.hpp"
#include "SharedPointer.hpp"
#include "impl/SharedPointer.hpp"
#include "tarch/logging/Log.h"
#include "utils/Helpers.hpp"
#include "tarch/la/DynamicColumnMatrix.h"
#include "boost/tuple/tuple.hpp" // ??

namespace precice { namespace cplscheme { namespace tests {
class SerialImplicitCouplingSchemeTest;
} } }

namespace precice {
namespace cplscheme {

/**
 * @brief Serial coupling scheme with iterations per timestep to achieve strong solution.
 */
class SerialImplicitCouplingScheme : public ImplicitCouplingScheme
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
   * @param maxIterations [IN] Maximal iterations per coupling timestep.
   * @param monitorIterations [IN] If true, a txt file monitoring iterations is
   *                          written.
   */
  SerialImplicitCouplingScheme (
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
  
  /**
   * @brief Advances within the coupling scheme (not necessarily in time).
   *
   * Preconditions:
   * - initialize() has been called.
   */
  virtual void advance();
  
  // @brief Logging device.
  static tarch::logging::Log _log;
  
  friend class tests::SerialImplicitCouplingSchemeTest;
};

}} // namespace precice, cplscheme

#endif /* PRECICE_CPLSCHEME_SERIALIMPLICITCOUPLINGSCHEME_HPP_ */
