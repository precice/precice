#pragma once

#include <string>
#include "BiCouplingScheme.hpp"
#include "cplscheme/BaseCouplingScheme.hpp"
#include "cplscheme/Constants.hpp"
#include "logging/Logger.hpp"
#include "m2n/SharedPointer.hpp"
#include "utils/assertion.hpp"

namespace precice {

namespace testing {
// Forward declaration to friend the boost test struct
struct ParallelCouplingSchemeFixture;
} // namespace testing

namespace cplscheme {

/**
 * @brief Coupling scheme for parallel coupling, i.e. simultaneous execution of two coupled participants
 *
 * For more information, look into Benjamin's thesis, Section 3.5.
 * https://mediatum.ub.tum.de/doc/1320661/document.pdf
 */
class ParallelCouplingScheme : public BiCouplingScheme {
  friend struct testing::ParallelCouplingSchemeFixture; // Make the fixture friend of this class
public:
  /**
   * @brief Constructor.
   *
   * @param[in] maxTime Simulation time limit, or UNDEFINED_MAX_TIME.
   * @param[in] maxTimeWindows Simulation time windows limit, or UNDEFINED_TIME_WINDOWS.
   * @param[in] timeWindowSize Simulation time window size.
   * @param[in] minTimeStepSize Minimum time step size.
   * @param[in] firstParticipant Name of participant starting simulation.
   * @param[in] secondParticipant Name of second participant in coupling.
   * @param[in] localParticipant Name of participant using this coupling scheme.
   * @param[in] m2n Communication object for com. between participants.
   * @param[in] dtMethod Method used for determining the time window size, see https://www.precice.org/couple-your-code-timestep-sizes.html
   * @param[in] cplMode Set implicit or explicit coupling
   * @param[in] maxIterations maximum number of coupling iterations allowed for implicit coupling per time window
   */
  ParallelCouplingScheme(
      double                        maxTime,
      int                           maxTimeWindows,
      double                        timeWindowSize,
      double                        minTimeStepSize,
      const std::string &           firstParticipant,
      const std::string &           secondParticipant,
      const std::string &           localParticipant,
      m2n::PtrM2N                   m2n,
      constants::TimesteppingMethod dtMethod,
      CouplingMode                  cplMode,
      int                           minIterations,
      int                           maxIterations);

  ParallelCouplingScheme(
      double                        maxTime,
      int                           maxTimeWindows,
      double                        timeWindowSize,
      double                        minTimeStepSize,
      const std::string &           firstParticipant,
      const std::string &           secondParticipant,
      const std::string &           localParticipant,
      m2n::PtrM2N                   m2n,
      constants::TimesteppingMethod dtMethod,
      CouplingMode                  cplMode);

private:
  logging::Logger _log{"cplscheme::ParallelCouplingScheme"};

  /// @copydoc cplscheme::BaseCouplingScheme::exchangeInitialData()
  void exchangeInitialData() override final;

  void exchangeFirstData() override final;

  void exchangeSecondData() override final;

  const DataMap &getAccelerationData() override final;
};

} // namespace cplscheme
} // namespace precice
