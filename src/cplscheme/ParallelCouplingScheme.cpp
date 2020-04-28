#include "ParallelCouplingScheme.hpp"
#include "acceleration/Acceleration.hpp"
#include "m2n/M2N.hpp"
#include "math/math.hpp"
#include "utils/EigenHelperFunctions.hpp"
#include "utils/MasterSlave.hpp"

namespace precice {
namespace cplscheme {

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
    int                           maxIterations)
    : BiCouplingScheme(maxTime, maxTimeWindows, timeWindowSize, validDigits, firstParticipant,
        secondParticipant, localParticipant, m2n, maxIterations, cplMode, dtMethod) {}

void ParallelCouplingScheme::checkConfiguration()
{
  if (isImplicitCouplingScheme()) {
    PRECICE_CHECK(not getSendData().empty(), "No send data configured. Use explicit scheme for one-way coupling.");
  }
}

void ParallelCouplingScheme::initializeImplementation()
{
  determineInitialSend(getSendData());
  determineInitialReceive(getReceiveData());
}

void ParallelCouplingScheme::exchangeInitialData()
{
  // F: send, receive, S: receive, send
  if (doesFirstStep()) {
    if (sendsInitializedData()) {
      sendData();
    }
    if (receivesInitializedData()) {
      receiveData();
    }
  } else { // second participant
    if (receivesInitializedData()) {
      receiveData();
      // second participant has to save values for extrapolation
      updateOldValues(getReceiveData());
    }
    if (sendsInitializedData()) {
      updateOldValues(getSendData());
      sendData();
    }
  }
}

bool ParallelCouplingScheme::exchangeDataAndAccelerate()
{
  bool convergence = true;

  if (doesFirstStep()) { //first participant
    PRECICE_DEBUG("Sending data...");
    sendData();
    PRECICE_DEBUG("Receiving data...");
    if(isImplicitCouplingScheme()) {
      convergence = receiveConvergence();
    }
    receiveData();
  } else { //second participant
    PRECICE_DEBUG("Receiving data...");
    receiveData();
    if (isImplicitCouplingScheme()) {
      PRECICE_DEBUG("Perform acceleration (only second participant)...");
      convergence = accelerate();
      sendConvergence(convergence);
    }
    PRECICE_DEBUG("Sending data...");
    sendData();
  }

  return convergence;
}

void ParallelCouplingScheme::mergeData()
{
  PRECICE_TRACE();
  PRECICE_ASSERT(!doesFirstStep(), "Only the second participant should do the acceleration.");
  PRECICE_ASSERT(_allData.empty(), "This function should only be called once.");
  _allData.insert(getSendData().begin(), getSendData().end());
  _allData.insert(getReceiveData().begin(), getReceiveData().end());
}

} // namespace cplscheme
} // namespace precice
