#pragma once

#include "BaseCouplingScheme.hpp"
#include "logging/Logger.hpp"

namespace precice {
namespace cplscheme {

/**
 * @brief TODO
 */
class MultiCouplingScheme : public BaseCouplingScheme {
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
 * @param[in] dtMethod method to determine timestep size, see https://github.com/precice/precice/wiki/Adapter's-Time-Step-Sizes
 * @param[in] maxIterations maximum number of coupling sub-iterations allowed.
 */
  MultiCouplingScheme(
      double                        maxTime,
      int                           maxTimeWindows,
      double                        timeWindowSize,
      int                           validDigits,
      const std::string &           localParticipant,
      std::vector<m2n::PtrM2N>      m2ns,
      constants::TimesteppingMethod dtMethod,
      int                           maxIterations = -1);

  /// Adds data to be sent on data exchange and possibly be modified during coupling iterations.
  void addDataToSend(
      mesh::PtrData data,
      mesh::PtrMesh mesh,
      bool          initialize,
      int           index);

  /// Adds data to be received on data exchange.
  void addDataToReceive(
      mesh::PtrData data,
      mesh::PtrMesh mesh,
      bool          initialize,
      int           index);

  /**
   * @brief TODO
   */
  void checkConfiguration() override;

private:
  CouplingData *getData(int dataID);

  /// Communication device to the other coupling participant.
  std::vector<m2n::PtrM2N> _m2ns;

  /// Map from data ID -> all data (receive and send) with that ID
  DataMap _allData;

  std::vector<DataMap> _receiveDataVector;
  std::vector<DataMap> _sendDataVector;

  logging::Logger _log{"cplscheme::MultiCouplingScheme"};

  /**
   * @brief TODO
   */
  bool exchangeDataAndAccelerate() override;

  DataMap &getAccelerationData() override
  {
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

  /**
   * @brief TODO
   */
  void assignDataToConvergenceMeasure(ConvergenceMeasure *convergenceMeasure, int dataID) override;

  /// @brief TODO
  void storeData() override {
    for (DataMap &sendData : _sendDataVector) {
      for (DataMap::value_type &pair : sendData) {
        if (pair.second->oldValues.size() > 0) {
          pair.second->oldValues.col(0) = *pair.second->values;
        }
      }
    }
    for (DataMap &receiveData : _receiveDataVector) {
      for (DataMap::value_type &pair : receiveData) {
        if (pair.second->oldValues.size() > 0) {
          pair.second->oldValues.col(0) = *pair.second->values;
        }
      }
    }

  }
};

} // namespace cplscheme
} // namespace precice
