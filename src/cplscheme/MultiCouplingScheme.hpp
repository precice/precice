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

namespace precice::cplscheme {
class CouplingData;
struct ExchangeData;

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
 * @param[in] maxTime Simulation time limit, or UNDEFINED_MAX_TIME.
 * @param[in] maxTimeWindows Simulation time windows limit, or UNDEFINED_TIME_WINDOWS.
 * @param[in] timeWindowSize Simulation time window size.
 * @param[in] localParticipant Name of participant using this coupling scheme.
 * @param[in] m2ns M2N communications to all other participants of coupling scheme.
 * @param[in] maxIterations maximum number of coupling sub-iterations allowed.
 */
  MultiCouplingScheme(
      double                             maxTime,
      int                                maxTimeWindows,
      double                             timeWindowSize,
      const std::string &                localParticipant,
      std::map<std::string, m2n::PtrM2N> m2ns,
      const std::string &                controller,
      int                                minIterations,
      int                                maxIterations);

  /// Adds data to be sent on data exchange and possibly be modified during coupling iterations.
  void addDataToSend(
      const mesh::PtrData &data,
      mesh::PtrMesh        mesh,
      bool                 requiresInitialization,
      bool                 exchangeSubsteps,
      const std::string &  to);

  /// Adds data to be received on data exchange.
  void addDataToReceive(
      const mesh::PtrData &data,
      mesh::PtrMesh        mesh,
      bool                 requiresInitialization,
      bool                 exchangeSubsteps,
      const std::string &  from);

  void determineInitialDataExchange() override;

  std::vector<std::string> getCouplingPartners() const final override;

  bool hasAnySendData() final override;

private:
  /**
   * @brief A vector of m2ns. A m2n is a communication device to the other coupling participant.
   */
  std::map<std::string, m2n::PtrM2N> _m2ns;

  /**
   * @brief A vector of all data to be received.
   */
  std::map<std::string, DataMap> _receiveDataVector;

  /**
   * @brief A vector of all data to be sent.
   */
  std::map<std::string, DataMap> _sendDataVector;

  /// Coupling partners to receive initial data from
  std::set<std::string> _receiveInitialFrom;

  /// Coupling partners to send initial data to
  std::set<std::string> _sendInitialTo;

  logging::Logger _log{"cplscheme::MultiCouplingScheme"};

  void exchangeFirstData() final override;

  void exchangeSecondData() final override;

  bool sendsInitializedDataTo(const std::string &to) const;

  bool receivesInitializedDataFrom(const std::string &from) const;

  DataMap &getAccelerationData() final override;

  /// @copydoc cplscheme::BaseCouplingScheme::initializeReceiveDataStorage()
  void initializeReceiveDataStorage() final override;

  /// @copydoc cplscheme::BaseCouplingScheme::exchangeInitialData()
  void exchangeInitialData() final override;

  /// name of the controller participant
  std::string _controller;

  /// if this is the controller or not
  bool _isController;
};

} // namespace precice::cplscheme
