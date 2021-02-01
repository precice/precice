#pragma once

#include <string>
#include <vector>
#include "BaseCouplingScheme.hpp"
#include "cplscheme/Constants.hpp"
#include "logging/Logger.hpp"
#include "m2n/SharedPointer.hpp"
#include "mesh/SharedPointer.hpp"

namespace precice {
namespace cplscheme {
struct CouplingData;

/**
 * @brief A coupling scheme with multiple participants.
 *
 * ! General description
 * A MultiCouplingScheme couples multiple participants in a fully implicit fashion.
 * It is a specialization of BaseCouplingScheme.
 *
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
 * @param[in] dtMethod Method used for determining the time window size, see https://www.precice.org/couple-your-code-timestep-sizes.html
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

  /// returns list of all coupling partners
  std::vector<std::string> getCouplingPartners() const override final;

  /**
   * @returns true, if coupling scheme has any sendData
   */
  bool hasAnySendData() override final
  {
    return std::any_of(_sendDataVector.cbegin(), _sendDataVector.cend(), [](DataMap sendData) { return not sendData.empty(); });
  }

private:
  /**
   * @brief get CouplingData from _allData using dataID
   * @param dataID identifies CouplingData to be searched for
   * @return CouplingData from _allData corresponding to dataID
   */
  CouplingData *getData(int dataID);

  /**
   * @brief A vector of m2ns. A m2n is a communication device to the other coupling participant.
   */
  std::vector<m2n::PtrM2N> _m2ns;

  /**
   * @brief Map from data ID -> all data (receive and send) with that ID
   */
  DataMap _allData;

  /**
   * @brief A vector of all data to be received.
   */
  std::vector<DataMap> _receiveDataVector;

  /**
   * @brief A vector of all data to be sent.
   */
  std::vector<DataMap> _sendDataVector;

  logging::Logger _log{"cplscheme::MultiCouplingScheme"};

  /**
   * @brief Exchanges all data between the participants of the MultiCouplingScheme and applies acceleration.
   * @returns true, if iteration converged
   */
  bool exchangeDataAndAccelerate() override;

  /**
   * @brief MultiCouplingScheme applies acceleration to _allData
   * @returns DataMap bein accelerated
   */
  DataMap &getAccelerationData() override
  {
    return _allData;
  }

  /**
   * @brief Initialization of MultiCouplingScheme is similar to ParallelCouplingScheme. We only have to iterate over all pieces of data in _sendDataVector and _receiveDataVector.
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

  /**
   * @brief Needed for setting up convergence measures
   * @param convMeasure Convergence measure to which the data field is assigned to
   * @param dataID Data field to be assigned
   */
  void assignDataToConvergenceMeasure(ConvergenceMeasureContext *convergenceMeasure, int dataID) override;

  /**
   * @brief MultiCouplingScheme has to call store for all receive and send data in the vectors
   */
  void storeData() override
  {
    for (DataMap &sendData : _sendDataVector) {
      store(sendData);
    }
    for (DataMap &receiveData : _receiveDataVector) {
      store(receiveData);
    }
  }
};

} // namespace cplscheme
} // namespace precice
