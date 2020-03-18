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

void SerialCouplingScheme::initializeImplicit()
{
  PRECICE_CHECK(not getSendData().empty(), "No send data configured! Use explicit scheme for one-way coupling.");
  if (not doesFirstStep()) {
    if (not getConvergenceMeasures().empty()) {
      setupConvergenceMeasures();              // needs _couplingData configured
      setupDataMatrices(getAcceleratedData()); // Reserve memory and initialize data with zero
    }
    if (getAcceleration().get() != nullptr) {
      getAcceleration()->initialize(getAcceleratedData()); // Reserve memory, initialize
    }
  } else if (getAcceleration().get() != nullptr && not getAcceleration()->getDataIDs().empty()) {
    int dataID = *(getAcceleration()->getDataIDs().begin());
    PRECICE_CHECK(getSendData(dataID) == nullptr,
                  "In case of serial coupling, acceleration can be defined for "
                      << "data of second participant only!");
  }
}

void SerialCouplingScheme::checkInitialize()
{
  if (anyDataRequiresInitialization(getSendData())) {
    PRECICE_CHECK(not doesFirstStep(), "Only second participant can initialize data and send it!");
  }

  if (anyDataRequiresInitialization(getReceiveData())) {
    PRECICE_CHECK(doesFirstStep(), "Only first participant can receive initial data!");
  }
}

void SerialCouplingScheme::initializeImplementation()
{
  checkInitialize();

  if (anyDataRequiresInitialization(getSendData())) {
    hasToSendInitializedData();
  }

  if (anyDataRequiresInitialization(getReceiveData())) {
    hasToReceiveInitializedData();
  }

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
      if (isImplicitCouplingScheme() && getExtrapolationOrder() > 0) {
        doExtrapolationOn(getSendData());
      }
      // The second participant sends the initialized data to the first participant
      // here, which receives the data on call of initialize().
      sendData(getM2N());
      receiveAndSetTimeWindowSize();
      // This receive replaces the receive in initialize().
      receiveData(getM2N());
    }
  }
}

std::pair<bool, bool> SerialCouplingScheme::doAdvance()
{
  bool convergence, convergenceCoarseOptimization, doOnlySolverEvaluation; // @todo having the bools for convergence measurement declared for explicit and implicit coupling is not nice

  // initialize advance
  if (isImplicitCouplingScheme()) {
    PRECICE_DEBUG("Computed full length of iteration");
    convergenceCoarseOptimization = true;
    doOnlySolverEvaluation        = false;
  }

  if (doesFirstStep()) { // first participant
    PRECICE_DEBUG("Sending data...");
    sendData(getM2N());
    if (isImplicitCouplingScheme()) {
      convergence = checkConvergence();
    }
    PRECICE_DEBUG("Receiving data...");
    if (isExplicitCouplingScheme()) {
      receiveAndSetTimeWindowSize();
    }
    receiveData(getM2N());
  } else { // second participant
    if (isImplicitCouplingScheme()) {
      PRECICE_DEBUG("Test Convergence and accelerate...");
      ValuesMap designSpecifications;  // TODO make this better?
      int       accelerationShift = 1; // TODO @BU: why do we need an "accelerationShift" for SerialCouplingScheme, but not for the ParallelCouplingScheme?
      doAcceleration(designSpecifications, convergence, convergenceCoarseOptimization, doOnlySolverEvaluation, accelerationShift);
      getM2N()->send(convergence);
      getM2N()->send(getIsCoarseModelOptimizationActive());
    }
    PRECICE_DEBUG("Sending data...");
    sendData(getM2N());
    // the second participant does not want new data in the last iteration of the last time window
    if (isCouplingOngoing() || (isImplicitCouplingScheme() && not convergence)) {
      PRECICE_DEBUG("Receiving data...");
      receiveAndSetTimeWindowSize();
      receiveData(getM2N());
    }
  }
  return std::pair<bool, bool>(convergence, convergenceCoarseOptimization);
}

} // namespace cplscheme
} // namespace precice
