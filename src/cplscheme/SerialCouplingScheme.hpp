#pragma once

#include <string>
#include "BaseCouplingScheme.hpp"
#include "BiCouplingScheme.hpp"
#include "cplscheme/Constants.hpp"
#include "logging/Logger.hpp"
#include "m2n/SharedPointer.hpp"

namespace precice {

namespace testing {
// Forward declaration to friend the boost test struct
struct SerialCouplingSchemeFixture;
} // namespace testing

namespace cplscheme {
/**
 * @brief Coupling scheme for serial coupling, i.e. staggered execution of two coupled participants
 *
 * For more information, look into Benjamin's thesis, Section 3.5.
 * https://mediatum.ub.tum.de/doc/1320661/document.pdf
 */
class SerialCouplingScheme : public BiCouplingScheme {
  friend struct testing::SerialCouplingSchemeFixture; // Make the fixture friend of this class
public:
  /**
 * @brief Constructor.
 *
 * @param[in] maxTime Simulation time limit, or UNDEFINED_MAX_TIME.
 * @param[in] maxTimeWindows Simulation time windows limit, or UNDEFINED_TIME_WINDOWS.
 * @param[in] timeWindowSize Simulation time window size.
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
      const std::string &           firstParticipant,
      const std::string &           secondParticipant,
      const std::string &           localParticipant,
      m2n::PtrM2N                   m2n,
      constants::TimesteppingMethod dtMethod,
      CouplingMode                  cplMode,
      int                           minIterations,
      int                           maxIterations);

  SerialCouplingScheme(
      double                        maxTime,
      int                           maxTimeWindows,
      double                        timeWindowSize,
      const std::string &           firstParticipant,
      const std::string &           secondParticipant,
      const std::string &           localParticipant,
      m2n::PtrM2N                   m2n,
      constants::TimesteppingMethod dtMethod,
      CouplingMode                  cplMode);

  ImplicitData implicitDataToReceive() const override final;

private:
  logging::Logger _log{"cplschemes::SerialCouplingSchemes"};

  /// Determines, if the time window size is set by the participant.
  bool _participantSetsTimeWindowSize = false;

  /// Determines, if the time window size is received by the participant.
  bool _participantReceivesTimeWindowSize = false;

  /// Sends time window size, if this participant is the one to send
  void sendTimeWindowSize();

  /// Receives and sets the time window size, if this participant is the one to receive
  void receiveAndSetTimeWindowSize();

  /// @copydoc cplscheme::BaseCouplingScheme::exchangeInitialData()
  void exchangeInitialData() override final;

  void exchangeFirstData() override final;

  void exchangeSecondData() override final;

  DataMap &getAccelerationData() override final;
};

} // namespace cplscheme
} // namespace precice
