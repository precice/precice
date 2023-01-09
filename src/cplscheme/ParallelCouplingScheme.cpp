#include "ParallelCouplingScheme.hpp"

#include <utility>

#include "cplscheme/BiCouplingScheme.hpp"
#include "logging/LogMacros.hpp"

namespace precice::cplscheme {

ParallelCouplingScheme::ParallelCouplingScheme(
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
    int                           maxIterations,
    int                           extrapolationOrder)
    : BiCouplingScheme(maxTime, maxTimeWindows, timeWindowSize, validDigits, firstParticipant,
                       secondParticipant, localParticipant, std::move(m2n), maxIterations, cplMode, dtMethod, extrapolationOrder) {}

void ParallelCouplingScheme::exchangeInitialData()
{
  // F: send, receive, S: receive, send

  if (doesFirstStep()) {
    if (sendsInitializedData()) {
      sendData(getM2N(), getSendData());
    }

    if (receivesInitializedData()) {
      receiveData(getM2N(), getReceiveData());
      checkDataHasBeenReceived();
    } else {
      initializeZeroReceiveData(getReceiveData());
    }
  } else { // second participant
    if (receivesInitializedData()) {
      receiveData(getM2N(), getReceiveData());
      checkDataHasBeenReceived();
    } else {
      initializeZeroReceiveData(getReceiveData());
    }

    if (sendsInitializedData()) {
      sendData(getM2N(), getSendData());
    }
  }
}

void ParallelCouplingScheme::exchangeFirstData()
{
  if (doesFirstStep()) { // first participant
    PRECICE_DEBUG("Sending data...");
    sendData(getM2N(), getSendData());
  } else { // second participant
    PRECICE_DEBUG("Receiving data...");
    receiveData(getM2N(), getReceiveData());
    checkDataHasBeenReceived();
  }
}

void ParallelCouplingScheme::exchangeSecondData()
{
  if (isImplicitCouplingScheme()) {
    if (doesFirstStep()) { // first participant
      receiveConvergence(getM2N());
    } else { // second participant
      PRECICE_DEBUG("Perform acceleration (only second participant)...");
      doImplicitStep();
      sendConvergence(getM2N());
    }
  }

  if (doesFirstStep()) {
    PRECICE_DEBUG("Receiving data...");
    receiveData(getM2N(), getReceiveData());
    checkDataHasBeenReceived();
  }

  if (hasConverged() || isExplicitCouplingScheme()) {
    for (const auto &data : _cplData | boost::adaptors::map_values) {
      data->moveTimeStepsStorage();
    }
  }
  if (isImplicitCouplingScheme()) {
    storeIteration();
  }

  if (!doesFirstStep()) { // second participant
    PRECICE_DEBUG("Sending data...");
    sendData(getM2N(), getSendData());
  }
}

const DataMap ParallelCouplingScheme::getAccelerationData()
{
  // ParallelCouplingScheme applies acceleration to all CouplingData
  PRECICE_ASSERT(!doesFirstStep(), "Only the second participant should do the acceleration.");
  return _cplData;
}

} // namespace precice::cplscheme
