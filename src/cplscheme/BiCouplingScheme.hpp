#pragma once

#include <string>
#include <vector>
#include "BaseCouplingScheme.hpp"
#include "cplscheme/Constants.hpp"
#include "logging/Logger.hpp"
#include "m2n/SharedPointer.hpp"
#include "mesh/SharedPointer.hpp"
#include "utils/assertion.hpp"

// Forward declaration to friend the boost test struct
namespace CplSchemeTests {
namespace SerialImplicitCouplingSchemeTests {
struct testExtrapolateData;
}
} // namespace CplSchemeTests

namespace precice {
namespace cplscheme {
struct CouplingData;

/**
 * @brief Abstract base class for coupling schemes with two participants.
 *
 * ! General description
 * A BiCouplingScheme couples two participants. It is a specialization of
 * BaseCouplingScheme.
 *
 */
class BiCouplingScheme : public BaseCouplingScheme {

public:
  BiCouplingScheme(
      double                        maxTime,
      int                           maxTimeWindows,
      double                        timeWindowSize,
      int                           validDigits,
      const std::string &           firstParticipant,
      const std::string &           secondParticipant,
      const std::string &           localParticipant,
      m2n::PtrM2N                   m2n,
      int                           maxIterations,
      CouplingMode                  cplMode,
      constants::TimesteppingMethod dtMethod);

  friend struct CplSchemeTests::SerialImplicitCouplingSchemeTests::testExtrapolateData; // For whitebox tests

  /// Adds data to be sent on data exchange and possibly be modified during coupling iterations.
  void addDataToSend(
      mesh::PtrData data,
      mesh::PtrMesh mesh,
      bool          requiresInitialization);

  /// Adds data to be received on data exchange.
  void addDataToReceive(
      mesh::PtrData data,
      mesh::PtrMesh mesh,
      bool          requiresInitialization);

  /// returns list of all coupling partners
  std::vector<std::string> getCouplingPartners() const override final;

  /**
   * @returns true, if coupling scheme has any sendData
   */
  bool hasAnySendData() override final
  {
    return not getSendData().empty();
  }

  /**
   * @returns true, if coupling scheme has sendData with given DataID
   */
  bool hasSendData(int dataID)
  {
    return getSendData(dataID) != nullptr;
  }

protected:
  /// Returns all data to be sent.
  DataMap &getSendData()
  {
    return _sendData;
  }

  /// Returns all data to be received.
  DataMap &getReceiveData()
  {
    return _receiveData;
  }

  /// Sets the values
  CouplingData *getSendData(int dataID);

  /// Returns all data to be received with data ID as given.
  CouplingData *getReceiveData(int dataID);

  /// @return Communication device to the other coupling participant.
  m2n::PtrM2N getM2N() const
  {
    PRECICE_ASSERT(_m2n);
    return _m2n;
  }

  /// @brief Receive from coupling partner and return whether coupling scheme has converged
  bool receiveConvergence();

private:
  mutable logging::Logger _log{"cplscheme::BiCouplingScheme"};

  /// Communication to the other coupling participant.
  m2n::PtrM2N _m2n;

  /// All send data as a map "data ID -> data"
  DataMap _sendData;

  /// All receive data as a map "data ID -> data"
  DataMap _receiveData;

  /// First participant name.
  std::string _firstParticipant = "unknown";

  /// Second participant name.
  std::string _secondParticipant = "unknown";

  /// Implements functionality for setupConvergenceMeasures
  void assignDataToConvergenceMeasure(ConvergenceMeasureContext *convMeasure, int dataID) override
  {
    if ((getSendData(dataID) != nullptr)) {
      convMeasure->couplingData = getSendData(dataID);
    } else {
      convMeasure->couplingData = getReceiveData(dataID);
      PRECICE_ASSERT(convMeasure->couplingData != nullptr);
    }
  }

  /**
   * @brief BiCouplingScheme has to call store for receive and send data
   */
  void storeData() override
  {
    store(getSendData());
    store(getReceiveData());
  }
};

} // namespace cplscheme
} // namespace precice
