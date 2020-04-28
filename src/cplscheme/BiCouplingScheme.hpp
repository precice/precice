#pragma once

#include "BaseCouplingScheme.hpp"
#include "logging/Logger.hpp"

namespace precice {
namespace cplscheme {

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

protected:
  /// Returns all data to be sent.
  DataMap &getSendData()
  {
    return _sendData;
  }

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

  /// TODO
  bool receiveConvergence();

private:

  mutable logging::Logger _log{"cplscheme::BiCouplingScheme"};

  /// Communication device to the other coupling participant.
  m2n::PtrM2N _m2n;

  /// Map from data ID -> all send data with that ID
  DataMap _sendData;

  /// Map from data ID -> all receive data with that ID
  DataMap _receiveData;

  /// Implements functionality for setupConvergenceMeasures
  void assignDataToConvergenceMeasure(ConvergenceMeasure* convMeasure, int dataID) override {
    if ((getSendData(dataID) != nullptr)) {
      convMeasure->couplingData = getSendData(dataID);
    } else {
      convMeasure->couplingData = getReceiveData(dataID);
      PRECICE_ASSERT(convMeasure->couplingData != nullptr);
    }
  }

  /// @brief TODO
  void storeData() override {
    for (DataMap::value_type &pair : getSendData()) {
      if (pair.second->oldValues.size() > 0) {
        pair.second->oldValues.col(0) = *pair.second->values;
      }
    }
    for (DataMap::value_type &pair : getReceiveData()) {
      if (pair.second->oldValues.size() > 0) {
        pair.second->oldValues.col(0) = *pair.second->values;
      }
    }
  }
};

} // namespace cplscheme
} // namespace precice
