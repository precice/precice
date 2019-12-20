#pragma once

#include "BaseCouplingScheme.hpp"
#include "logging/Logger.hpp"

// Forward declaration to friend the boost test struct
namespace CplSchemeTests {
namespace SerialImplicitCouplingSchemeTests {
struct testExtrapolateData;
}
} // namespace CplSchemeTests

namespace precice {
namespace cplscheme {

/**
 * @brief Coupling scheme for serial coupling, i.e. staggered execution of two coupled participants
 *
 * For more information, look into Benjamin's thesis, Section 3.5. 
 * https://mediatum.ub.tum.de/doc/1320661/document.pdf
 */
class SerialCouplingScheme : public BaseCouplingScheme {
public:
  /**
 * @brief Constructor.
 *
 * @param[in] maxTime Simulation time limit, or UNDEFINED_TIME.
 * @param[in] maxTimesteps Simulation timestep limit, or UNDEFINED_TIMESTEPS.
 * @param[in] timestepLength Simulation timestep length.
 * @param[in] firstParticipant Name of participant starting simulation.
 * @param[in] secondParticipant Name of second participant in coupling.
 * @param[in] localParticipant Name of participant using this coupling scheme.
 * @param[in] communication Communication object for com. between participants.
 * @param[in] monitorIterations If true, a txt file monitoring iterations is written.
 */
  SerialCouplingScheme(
      double                        maxTime,
      int                           maxTimesteps,
      double                        timestepLength,
      int                           validDigits,
      const std::string &           firstParticipant,
      const std::string &           secondParticipant,
      const std::string &           localParticipant,
      m2n::PtrM2N                   m2n,
      constants::TimesteppingMethod dtMethod,
      CouplingMode                  cplMode,
      int                           maxIterations = 1);

  virtual void initialize(double startTime, int startTimestep);

  virtual void initializeData();

  virtual void advance();

  logging::Logger _log{"cplschemes::SerialCouplingSchemes"};

  friend struct CplSchemeTests::SerialImplicitCouplingSchemeTests::testExtrapolateData; // For whitebox tests
};

} // namespace cplscheme
} // namespace precice
