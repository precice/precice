#include "SerialCouplingScheme.hpp"
#include "acceleration/Acceleration.hpp"
#include "m2n/M2N.hpp"
#include "math/math.hpp"
#include "utils/EigenHelperFunctions.hpp"
#include "utils/MasterSlave.hpp"

namespace precice {
namespace cplscheme {

SerialCouplingScheme::SerialCouplingScheme(
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
    : BaseCouplingScheme(maxTime, maxTimeWindows, timeWindowSize, validDigits, firstParticipant,
                         secondParticipant, localParticipant, m2n, maxIterations, cplMode, dtMethod) {}

void SerialCouplingScheme::checkForSend()
{
  PRECICE_CHECK(not getSendData().empty(), "No send data configured! Use explicit scheme for one-way coupling.");
}

void SerialCouplingScheme::checkAcceleration()
{
  if (doesFirstStep() && getAcceleration() && not getAcceleration()->getDataIDs().empty()) {
    int dataID = *(getAcceleration()->getDataIDs().begin());
    PRECICE_CHECK(getSendData(dataID) == nullptr,
                  "In case of serial coupling, acceleration can be defined for "
                      << "data of second participant only!");
  }
}

void SerialCouplingScheme::initializeImplementation()
{
  // Perform checks for initialization of serial coupling.
  if (anyDataRequiresInitialization(getSendData())) {
    PRECICE_CHECK(not doesFirstStep(), "In serial coupling only second participant can initialize data and send it!");
  }

  if (anyDataRequiresInitialization(getReceiveData())) {
    PRECICE_CHECK(doesFirstStep(), "In serial coupling only first participant can receive initial data!");
  }

  // Initialize coupling
  determineInitialSend(getSendData());
  determineInitialReceive(getReceiveData());

  // If the second participant initializes data, the first receive for the
  // second participant is done in initializeData() instead of initialize().
  if (not doesFirstStep() && not sendsInitializedData() && isCouplingOngoing()) {
    PRECICE_DEBUG("Receiving data");
    receiveAndSetTimeWindowSize();
    receiveData(getM2N());
  }
}

void SerialCouplingScheme::exchangeInitialData()
{
  PRECICE_ASSERT(isCouplingOngoing());
  if (doesFirstStep()) {
    PRECICE_ASSERT(not sendsInitializedData(), "First participant cannot send data during initialization.");
    if (receivesInitializedData()) {
      receiveData(getM2N());
    }
  } else { // second participant
    PRECICE_ASSERT(not receivesInitializedData(), "Only first participant can receive data during initialization.");
    if (sendsInitializedData()) {
      updateOldValues(getSendData());
      // The second participant sends the initialized data to the first participant
      // here, which receives the data on call of initialize().
      sendData(getM2N());
      receiveAndSetTimeWindowSize();
      // This receive replaces the receive in initialize().
      receiveData(getM2N());
    }
  }
}

std::pair<bool, bool> SerialCouplingScheme::exchangeDataAndAccelerate()
{
  bool convergence, convergenceCoarseOptimization; // @todo having the bools for convergence measurement declared for explicit and implicit coupling is not nice

  if (doesFirstStep()) { // first participant
    PRECICE_DEBUG("Sending data...");
    sendData(getM2N());
    if (isImplicitCouplingScheme()) {
      convergence = receiveConvergence();
      if (convergence) {
        timeWindowCompleted();
      }
    }
    PRECICE_DEBUG("Receiving data...");
    receiveData(getM2N());
  } else { // second participant
    if (isImplicitCouplingScheme()) {
      PRECICE_DEBUG("Test Convergence and accelerate...");
      int       accelerationShift = 1; // TODO @BU: why do we need an "accelerationShift" for SerialCouplingScheme, but not for the ParallelCouplingScheme?
      std::pair<bool, bool> convergenceInformation = accelerate(accelerationShift);
      convergence = convergenceInformation.first;
      convergenceCoarseOptimization = convergenceInformation.second;
    }
    PRECICE_DEBUG("Sending data...");
    sendData(getM2N());
    // the second participant does not want new data in the last iteration of the last time window
    if (isCouplingOngoing() || (isImplicitCouplingScheme() && not convergence)) {
      PRECICE_DEBUG("Receiving data...");
      receiveData(getM2N());
    }
  }
  return std::pair<bool, bool>(convergence, convergenceCoarseOptimization);
}

} // namespace cplscheme
} // namespace precice
