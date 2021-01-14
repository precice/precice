#pragma once

#include <string>
#include "BiCouplingScheme.hpp"
#include "cplscheme/BaseCouplingScheme.hpp"
#include "cplscheme/Constants.hpp"
#include "logging/Logger.hpp"
#include "m2n/SharedPointer.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace cplscheme {

/**
 * @brief Coupling scheme for parallel coupling, i.e. simultaneous execution of two coupled participants
 *
 * For more information, look into Benjamin's thesis, Section 3.5. 
 * https://mediatum.ub.tum.de/doc/1320661/document.pdf
 */
class ParallelCouplingScheme : public BiCouplingScheme {
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
  ParallelCouplingScheme(
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
  logging::Logger _log{"cplscheme::ParallelCouplingScheme"};

  /// Map from data ID -> all data (receive and send) with that ID
  DataMap _allData;

  /**
   * @brief Exchanges all data between the participants of the ParallelCouplingScheme and applies acceleration.
   * @returns true, if iteration converged
   */
  bool exchangeDataAndAccelerate() override;

  /**
   * @brief ParallelCouplingScheme applies acceleration to _allData
   * @returns DataMap being accelerated
   */
  DataMap &getAccelerationData() override
  {
    PRECICE_ASSERT(!doesFirstStep(), "Only the second participant should do the acceleration.");
    return _allData;
  }

  /**
   * @brief determine whether data has to be sent/received
   */
  void initializeImplementation() override;

  /**
   * @brief merges send and receive data into one map (for parallel acceleration)
   */
  void mergeData() override;

  /**
   * @brief Exchanges data, if it has to be initialized.
   */
  void exchangeInitialData() override;
};

} // namespace cplscheme
} // namespace precice
