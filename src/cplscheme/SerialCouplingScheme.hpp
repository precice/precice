#pragma once

#include <string>
#include "BaseCouplingScheme.hpp"
#include "BiCouplingScheme.hpp"
#include "cplscheme/Constants.hpp"
#include "logging/Logger.hpp"
#include "m2n/SharedPointer.hpp"

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
 * @param[in] validDigits valid digits for computation of the remainder of a time window
 * @param[in] firstParticipant Name of participant starting simulation.
 * @param[in] secondParticipant Name of second participant in coupling.
 * @param[in] localParticipant Name of participant using this coupling scheme.
 * @param[in] m2n Communication object for com. between participants.
 * @param[in] dtMethod Method used for determining the time window size, see https://www.precice.org/couple-your-code-timestep-sizes.html
 * @param[in] cplMode Set implicit or explicit coupling
 * @param[in] maxIterations maximum number of coupling iterations allowed for implicit coupling per time window
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

private:
  logging::Logger _log{"cplschemes::SerialCouplingSchemes"};

  /// Determines, if the time window size is set by the participant.
  bool _participantSetsTimeWindowSize = false;

  /// Determines, if the time window size is received by the participant.
  bool _participantReceivesTimeWindowSize = false;

  /// Receives and sets the time window size, if this participant is the one to receive
  void receiveAndSetTimeWindowSize();

  /**
   * @brief Exchanges data between the participants of the SerialCouplingSchemes and applies acceleration.
   * @returns true, if iteration converged
   */
  bool exchangeDataAndAccelerate() override;

  /**
   * @brief SerialCouplingSchemes applies acceleration to send data
   * @returns DataMap being accelerated
   */
  DataMap &getAccelerationData() override
  {
    return getSendData();
  }

  /**
   * @brief determine whether data has to be sent/received
   */
  void initializeImplementation() override;

  /**
   * @brief noop for SerialCouplingScheme
   */
  void mergeData() override{};

  /**
   * @brief Exchanges data, if it has to be initialized.
   */
  void exchangeInitialData() override;
};

} // namespace cplscheme
} // namespace precice
