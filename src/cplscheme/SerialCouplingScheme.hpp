#pragma once

#include "BaseCouplingScheme.hpp"
#include "BiCouplingScheme.hpp"
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
class SerialCouplingScheme : public BiCouplingScheme {
public:
  /**
 * @brief Constructor.
 *
 * @param[in] maxTime Simulation time limit, or UNDEFINED_TIME.
 * @param[in] maxTimeWindows Simulation time windows limit, or UNDEFINED_TIMEWINDOWS.
 * @param[in] timeWindowSize Simulation time window size.
 * @param[in] validDigits TODO
 * @param[in] firstParticipant Name of participant starting simulation.
 * @param[in] secondParticipant Name of second participant in coupling.
 * @param[in] localParticipant Name of participant using this coupling scheme.
 * @param[in] m2n Communication object for com. between participants. TODO?
 * @param[in] dtMethod TODO
 * @param[in] cplMode TODO
 * @param[in] maxIterations TODO
 */
  SerialCouplingScheme(
      double                        maxTime,
      int                           maxTimeWindows,
      double                        timeWindowSize,
      int                           validDigits,
      const std::string &           firstParticipant,
      const std::string &           secondParticipant,
      const std::string &           localParticipant,
      m2n::PtrM2N                   m2n,
      constants::TimesteppingMethod dtMethod,
      CouplingMode                  cplMode,
      int                           maxIterations = -1);

  friend struct CplSchemeTests::SerialImplicitCouplingSchemeTests::testExtrapolateData; // For whitebox tests

  /**
   * @brief TODO
   */
  void checkConfiguration() override;

private:

  /**
   * @brief TODO
   */
  logging::Logger _log{"cplschemes::SerialCouplingSchemes"};

  /// Determines, if the time window size is set by the participant.
  bool _participantSetsTimeWindowSize = false;

  /// Determines, if the time window size is received by the participant.
  bool _participantReceivesTimeWindowSize = false;

  /// Receives and sets the time window size, if this participant is the one to receive
  void receiveAndSetTimeWindowSize();

  /**
   * @brief TODO
   */
  bool exchangeDataAndAccelerate() override;

  DataMap &getAccelerationData() override
  {
    return getSendData();
  }

  /**
   * @brief TODO
   */
  void initializeImplementation() override;

  /**
   * @brief noop for SerialCouplingScheme
   */
  void mergeData() override {};

  /**
   * @brief TODO
   */
  void exchangeInitialData() override;
};

} // namespace cplscheme
} // namespace precice
