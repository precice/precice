// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_CPLSCHEME_PARALLELEXPLICITCOUPLINGSCHEME_HPP_
#define PRECICE_CPLSCHEME_PARALLELEXPLICITCOUPLINGSCHEME_HPP_

#include "ExplicitCouplingScheme.hpp"
#include "BaseCouplingScheme.hpp"
#include "Constants.hpp"
#include "io/TXTTableWriter.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/PropertyContainer.hpp"
#include "mesh/Mesh.hpp"
#include "tarch/logging/Log.h"
#include "utils/Helpers.hpp"
#include "boost/tuple/tuple.hpp"

namespace precice { namespace cplscheme { namespace tests {
//class ParallelExplicitCouplingSchemeTest;
} } }

namespace precice {
namespace cplscheme {

/// @brief Serial coupling scheme without iterations per timestep.
class ParallelExplicitCouplingScheme : public ExplicitCouplingScheme
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
  ParallelExplicitCouplingScheme (
    double                maxTime,
    int                   maxTimesteps,
    double                timestepLength,
    int                   validDigits,
    const std::string&    firstParticipant,
    const std::string&    secondParticipant,
    const std::string&    localParticipant,
    com::PtrCommunication communication,
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
   * @brief Advances within the coupling scheme.
   *
   * Preconditions:
   * - initialize() has been called.
   */
  virtual void advance();

  /// @brief Logging device.
  static tarch::logging::Log _log;

  //friend class tests::ParallelExplicitCouplingSchemeTest;
};

}} // namespace precice, cplscheme

#endif /* PRECICE_CPLSCHEME_PARALLELEXPLICITCOUPLINGSCHEME_HPP_ */
