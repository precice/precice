#pragma once

#include <map>
#include <string>
#include <vector>
#include "BaseCouplingScheme.hpp"
#include "cplscheme/Constants.hpp"
#include "logging/Logger.hpp"
#include "m2n/SharedPointer.hpp"
#include "mesh/SharedPointer.hpp"
#include "utils/Helpers.hpp"

namespace precice {

namespace testing {
// Forward declaration to friend the boost test struct
struct ParallelCouplingSchemeFixture;
} // namespace testing

namespace cplscheme {
class CouplingData;
struct ExchangeData;

/**
 * @brief A coupling scheme with multiple participants.
 *
 * ! General description
 * A ParallelCouplingScheme couples multiple participants in a fully implicit fashion.
 * It is a specialization of BaseCouplingScheme.
 *
 */
class ParallelCouplingScheme : public BaseCouplingScheme {
  friend struct testing::ParallelCouplingSchemeFixture; // Make the fixture friend of this class
public:
  /**
 * @brief Constructor.
 *
 * @param[in] maxTime Simulation time limit, or UNDEFINED_TIME.
 * @param[in] maxTimeWindows Simulation time windows limit, or UNDEFINED_TIMEWINDOWS.
 * @param[in] timeWindowSize Simulation time window size.
 * @param[in] validDigits valid digits for computation of the remainder of a time window
 * @param[in] localParticipant Name of participant using this coupling scheme.
 * @param[in] m2ns M2N communications to all other participants of coupling scheme.
 * @param[in] dtMethod Method used for determining the time window size, see https://www.precice.org/couple-your-code-timestep-sizes.html
 * @param[in] cplMode Set implicit or explicit coupling
 * @param[in] maxIterations maximum number of coupling sub-iterations allowed.
 * @param[in] extrapolationOrder order used for extrapolation
 */
  ParallelCouplingScheme(
      double                             maxTime,
      int                                maxTimeWindows,
      double                             timeWindowSize,
      int                                validDigits,
      const std::string &                localParticipant,
      std::map<std::string, m2n::PtrM2N> m2ns,
      constants::TimesteppingMethod      dtMethod,
      CouplingMode                       cplMode,
      const std::string &                controller,
      int                                maxIterations      = UNDEFINED_MAX_ITERATIONS,
      int                                extrapolationOrder = UNDEFINED_MAX_ITERATIONS);

private:
  logging::Logger _log{"cplscheme::ParallelCouplingScheme"};

  /**
   * @brief Exchanges all data between the participants of the ParallelCouplingScheme and applies acceleration.
   * @returns true, if iteration converged
   */
  bool exchangeDataAndAccelerate() override;

  /**
   * @brief ParallelCouplingScheme applies acceleration to all CouplingData
   * @returns DataMap being accelerated
   */
  const DataMap getAccelerationData() override;

  /**
   * @brief Exchanges data, if it has to be initialized.
   */
  void exchangeInitialData() override;
};

} // namespace cplscheme
} // namespace precice
