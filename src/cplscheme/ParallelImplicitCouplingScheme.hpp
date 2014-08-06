// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_CPLSCHEME_PARALLELIMPLICITCOUPLINGSCHEME_HPP_
#define PRECICE_CPLSCHEME_PARALLELIMPLICITCOUPLINGSCHEME_HPP_

#include "ParallelCouplingScheme.hpp"
#include "Constants.hpp"
#include "SharedPointer.hpp"
#include "tarch/logging/Log.h"

namespace precice { namespace cplscheme { namespace tests {
class ParallelImplicitCouplingSchemeTest;
} } }

namespace precice {
namespace cplscheme {

/// @brief Parallel coupling scheme which lets the participants run in parallel to each other.
class ParallelImplicitCouplingScheme : public ParallelCouplingScheme
{
public:
/**
 *
 * @brief Constructor.
 *
 * @param maxTime [IN] Simulation time limit, or UNDEFINED_TIME.
 * @param maxTimesteps [IN] Simulation timestep limit, or UNDEFINED_TIMESTEPS.
 * @param timestepLength [IN] Simulation timestep length.
 * @param firstParticipant [IN] Name of first participant.
 * @param secondParticipant [IN] Name of second participant who does the pp.
 * @param localParticipant [IN] Name of participant using this coupling scheme.
 * @param communication [IN] Communication object for com. between participants.
 * @param maxIterations [IN] Maximal iterations per coupling timestep.
 * @param monitorIterations [IN] If true, a txt file monitoring iterations is
 *                          written.
 */
  // call superconstructor in implementation
  ParallelImplicitCouplingScheme (
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

  /// @brief Initializes the coupling scheme.
  virtual void initialize (
    double startTime,
    int    startTimestep );

  /**
   * @brief Initializes data with written values.
   *
   * Preconditions:
   * - initialize() has been called.
   * - advance() has NOT yet been called.
   */
  virtual void initializeData();

  /**
   * @brief Advances within the coupling scheme (not necessarily in time).
   *
   * Preconditions:
   * - initialize() has been called.
   */
  virtual void advance();

  virtual std::string printCouplingState() const;
  
private:

  /// @brief Logging device.
  static tarch::logging::Log _log;

  /// @brief Map from data ID -> all data (receive and send) with that ID
  DataMap _allData;

  friend class tests::ParallelImplicitCouplingSchemeTest;
};

}} // namespace precice, cplscheme

#endif /* PRECICE_CPLSCHEME_PARALLELIMPLICITCOUPLINGSCHEME_HPP_ */
