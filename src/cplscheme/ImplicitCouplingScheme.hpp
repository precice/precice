// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_CPLSCHEME_IMPLICITCOUPLINGSCHEME_HPP_
#define PRECICE_CPLSCHEME_IMPLICITCOUPLINGSCHEME_HPP_

#include "BaseCouplingScheme.hpp"
#include "SharedPointer.hpp"
#include "Constants.hpp"
#include "impl/SharedPointer.hpp"
#include "io/TXTTableWriter.hpp"
#include "tarch/logging/Log.h"
#include "utils/Helpers.hpp"
#include "tarch/la/DynamicColumnMatrix.h"
#include "boost/tuple/tuple.hpp"

namespace precice { namespace cplscheme { namespace tests {
class ImplicitCouplingSchemeTest;
} } }

namespace precice {
namespace cplscheme {

/**
 * Abstract class that provides the basic functionalities for implicit coupling,
 * i.e. subiterating in every timestep to converge towards the strong solution.
 * The functionalities that differ for the classcial serial coupling and the
 * parallel coupling are implemented in the subclasses SerialImplicitCouplingScheme
 * and ParallelImplicitCouplingScheme.
 * Please look at ./coupling_steering.pdf for a brief sketch of the differences
 * between the serial and parallel implicit coupling.
 *
 * @brief Abstract coupling scheme with iterations per timestep to achieve strong solution.
 */
class ImplicitCouplingScheme : public BaseCouplingScheme
{
public:
  
  /**
   * @brief Constructor.
   *
   *
   * @param maxTime [IN] Simulation time limit, or UNDEFINED_TIME.
   * @param maxTimesteps [IN] Simulation timestep limit, or UNDEFINED_TIMESTEPS.
   * @param timestepLength [IN] Simulation timestep length.
   * @param firstParticipant [IN] Name of first participant in coupling.
   * @param secondParticipant [IN] Name of second participant in coupling.
   * @param localParticipant [IN] Name of participant using this coupling scheme.
   * @param communication [IN] Communication object for com. between participants.
   * @param maxIterations [IN] Maximal iterations per coupling timestep.
   * @param monitorIterations [IN] If true, a txt file monitoring iterations is
   *                          written.
   */
  ImplicitCouplingScheme (
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
  
 
  
//  virtual std::string printCouplingState() const;


private:

  // typedef tarch::la::DynamicColumnMatrix<double> DataMatrix;
  
  // typedef tarch::la::DynamicVector<double> DataVector;
  
  /// @brief Logging device.
  static tarch::logging::Log _log;
  
// friend class tests::ImplicitCouplingSchemeTest;
};

}} // namespace precice, cplscheme

#endif /* PRECICE_CPLSCHEME_IMPLICITCOUPLINGSCHEME_HPP_ */
