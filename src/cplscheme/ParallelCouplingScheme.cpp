#include "ParallelCouplingScheme.hpp"
#include <algorithm>
#include <cstddef>
#include <map>
#include <memory>
#include <ostream>
#include <type_traits>
#include <utility>
#include "acceleration/Acceleration.hpp"
#include "acceleration/SharedPointer.hpp"
#include "cplscheme/BaseCouplingScheme.hpp"
#include "cplscheme/CouplingData.hpp"
#include "cplscheme/SharedPointer.hpp"
#include "logging/LogMacros.hpp"
#include "m2n/SharedPointer.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"

namespace precice::cplscheme {

ParallelCouplingScheme::ParallelCouplingScheme(
    double                             maxTime,
    int                                maxTimeWindows,
    double                             timeWindowSize,
    int                                validDigits,
    const std::string &                localParticipant,
    std::map<std::string, m2n::PtrM2N> m2ns,
    constants::TimesteppingMethod      dtMethod,
    CouplingMode                       cplMode,
    const std::string &                controller,
    int                                maxIterations,
    int                                extrapolationOrder)
    : BaseCouplingScheme(maxTime, maxTimeWindows, timeWindowSize, validDigits, localParticipant, controller, maxIterations, cplMode, dtMethod, extrapolationOrder)
{
  _m2ns = m2ns;
  if (_m2ns.size() > 1) {
    // @todo implement ParallelCouplingScheme with multi coupling for explicit coupling
    PRECICE_ASSERT(isImplicitCouplingScheme(), "ParallelCouplingScheme with multi coupling is always Implicit.");
  }
  PRECICE_DEBUG("MultiCoupling scheme is created for {}.", localParticipant);
}

void ParallelCouplingScheme::exchangeInitialData()
{
  if (!doesFirstStep()) {
    if (receivesInitializedData()) {
      for (auto &receiveExchange : _receiveDataVector) {
        receiveData(_m2ns[receiveExchange.first], receiveExchange.second);
      }
      checkDataHasBeenReceived();
    }
    if (sendsInitializedData()) {
      for (auto &sendExchange : _sendDataVector) {
        sendData(_m2ns[sendExchange.first], sendExchange.second);
      }
    }
  } else {
    if (sendsInitializedData()) {
      for (auto &sendExchange : _sendDataVector) {
        sendData(_m2ns[sendExchange.first], sendExchange.second);
      }
    }
    if (receivesInitializedData()) {
      for (auto &receiveExchange : _receiveDataVector) {
        receiveData(_m2ns[receiveExchange.first], receiveExchange.second);
      }
      checkDataHasBeenReceived();
    }
  }
  PRECICE_DEBUG("Initial data is exchanged in ParallelCouplingScheme");
}

bool ParallelCouplingScheme::exchangeDataAndAccelerate()
{
  PRECICE_DEBUG("Computed full length of iteration");

  bool convergence = true;

  if (!doesFirstStep()) {
    PRECICE_DEBUG("Receiving data...");
    for (auto &receiveExchange : _receiveDataVector) {
      receiveData(_m2ns[receiveExchange.first], receiveExchange.second);
    }
    checkDataHasBeenReceived();
    if (isImplicitCouplingScheme()) {
      PRECICE_DEBUG("Perform acceleration (only second participant)...");
      convergence = doImplicitStep();
      for (const auto &m2nPair : _m2ns) {
        sendConvergence(m2nPair.second, convergence);
      }
    }
    PRECICE_DEBUG("Sending data...");
    for (auto &sendExchange : _sendDataVector) {
      sendData(_m2ns[sendExchange.first], sendExchange.second);
    }
  } else {
    PRECICE_DEBUG("Sending data...");
    for (auto &sendExchange : _sendDataVector) {
      sendData(_m2ns[sendExchange.first], sendExchange.second);
    }
    PRECICE_DEBUG("Receiving data...");
    if (isImplicitCouplingScheme()) {
      convergence = receiveConvergence(_m2ns[_controller]);
    }
    for (auto &receiveExchange : _receiveDataVector) {
      receiveData(_m2ns[receiveExchange.first], receiveExchange.second);
    }
    checkDataHasBeenReceived();
  }

  return convergence;
}

typedef std::map<int, PtrCouplingData> DataMap;

const DataMap ParallelCouplingScheme::getAccelerationData()
{
  DataMap accelerationData;
  for (auto &data : allCouplingData()) {
    PRECICE_ASSERT(accelerationData.count(data->getDataID()) == 0);
    accelerationData[data->getDataID()] = data;
  }
  return accelerationData;
}

} // namespace precice::cplscheme
