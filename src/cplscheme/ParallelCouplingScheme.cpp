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
  bool initialCommunication = true;

  if (doesFirstStep()) {
    if (sendsInitializedData()) {
      sendData(getM2N(), getSendData(), initialCommunication);
    } else {
      // need to clear storage, because otherwise send data will be polluted in next iteration
      for (const auto &data : getSendData() | boost::adaptors::map_values) {
        bool keepWindowStart = false;
        data->clearTimeStepsStorage(keepWindowStart);
      }
    }
    if (receivesInitializedData()) {
      receiveData(getM2N(), getReceiveData(), initialCommunication);
      storeReceiveData(time::Storage::WINDOW_END); // use constant data in first iteration
      checkDataHasBeenReceived();
    } else {
      initializeZeroReceiveData(getReceiveData());
    }
  } else { // second participant
    if (receivesInitializedData()) {
      receiveData(getM2N(), getReceiveData(), initialCommunication);
      storeReceiveData(time::Storage::WINDOW_END); // use constant data in first iteration
      checkDataHasBeenReceived();
    } else {
      initializeZeroReceiveData(getReceiveData());
    }
    if (sendsInitializedData()) {
      sendData(getM2N(), getSendData(), initialCommunication);
    } else {
      // need to clear storage, because otherwise send data will be polluted in next iteration
      for (const auto &data : getSendData() | boost::adaptors::map_values) {
        bool keepWindowStart = false;
        data->clearTimeStepsStorage(keepWindowStart);
      }
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
  if (doesFirstStep()) { // first participant
    PRECICE_DEBUG("Receiving data...");
    if (isImplicitCouplingScheme()) {
      receiveConvergence(getM2N());
    }
    receiveData(getM2N(), getReceiveData());
    checkDataHasBeenReceived();
  } else { // second participant
    if (isImplicitCouplingScheme()) {
      PRECICE_DEBUG("Perform acceleration (only second participant)...");
      doImplicitStep();
      sendConvergence(getM2N());
    }
    PRECICE_DEBUG("Sending data...");
    sendData(getM2N(), getSendData());
  }
  if (hasConverged()) {
    // @todo similar code breaks in SerialCouplingScheme.cpp for CplSchemeTests/SerialImplicitCouplingSchemeTests/ testConfiguredAbsConvergenceMeasureSynchronized. Why? @fsimonis
    // for (const auto &data : getAllData() | boost::adaptors::map_values) {
    //   data->moveTimeStepsStorage();
    // }
    for (const DataMap::value_type &pair : getReceiveData()) { // @todo this is probably an error which causes wrong initial data at the beginning of the window. But with getAllData() it also breaks...
      pair.second->moveTimeStepsStorage();
    }
  }
  if (isImplicitCouplingScheme()) {
    storeIteration();
  }
}

const DataMap ParallelCouplingScheme::getAccelerationData()
{
  // ParallelCouplingScheme applies acceleration to all CouplingData
  PRECICE_ASSERT(!doesFirstStep(), "Only the second participant should do the acceleration.");
  return getAllData();
}

} // namespace precice::cplscheme
