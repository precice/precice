#pragma once

#include "BiCouplingScheme.hpp"
#include "logging/Logger.hpp"

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
 * @param[in] validDigits TODO
 * @param[in] firstParticipant Name of participant starting simulation.
 * @param[in] secondParticipant Name of second participant in coupling.
 * @param[in] localParticipant Name of participant using this coupling scheme.
 * @param[in] m2n Communication object for com. between participants. TODO?
 * TODO add dtMethod, cplMode, maxIterations
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

  /**
   * @brief TODO
   */
  void checkConfiguration() override;

private:
  logging::Logger _log{"cplscheme::ParallelCouplingScheme"};

  /// Map from data ID -> all data (receive and send) with that ID
  DataMap _allData;

  /**
   * @brief TODO
   */
  bool exchangeDataAndAccelerate() override;

  DataMap &getAccelerationData() override
  {
    PRECICE_ASSERT(!doesFirstStep(), "Only the second participant should do the acceleration.");
    return _allData;
  }

  /**
   * @brief TODO
   */
  void initializeImplementation() override;

  /**
   * @brief merges send and receive data into one map (for parallel acceleration)
   */
  void mergeData() override;

  /**
   * @brief TODO
   */
  void exchangeInitialData() override;
};

} // namespace cplscheme
} // namespace precice
