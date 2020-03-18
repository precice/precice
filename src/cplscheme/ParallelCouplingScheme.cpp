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
    double                        timeWindowsSize,
    int                           validDigits,
    const std::string &           firstParticipant,
    const std::string &           secondParticipant,
    const std::string &           localParticipant,
    m2n::PtrM2N                   m2n,
    constants::TimesteppingMethod dtMethod,
    CouplingMode                  cplMode,
    int                           maxIterations)
    : BaseCouplingScheme(maxTime, maxTimeWindows, timeWindowsSize, validDigits, firstParticipant,
                         secondParticipant, localParticipant, m2n, maxIterations, cplMode, dtMethod) {}

void ParallelCouplingScheme::initializeImplicit()
{
  PRECICE_CHECK(not getSendData().empty(), "No send data configured. Use explicit scheme for one-way coupling.");
  if (not doesFirstStep()) {                 // second participant
    setupConvergenceMeasures();              // needs _couplingData configured
    mergeData();                             // merge send and receive data for all pp calls
    setupDataMatrices(getAcceleratedData()); // Reserve memory and initialize data with zero
    if (getAcceleration().get() != nullptr) {
      getAcceleration()->initialize(getAcceleratedData()); // Reserve memory, initialize
    }
  }
}

void ParallelCouplingScheme::initializeImplementation()
{
  if (anyDataRequiresInitialization(getSendData())) {
    hasToSendInitializedData();
  }

  if (anyDataRequiresInitialization(getReceiveData())) {
    hasToReceiveInitializedData();
  }
}

void ParallelCouplingScheme::exchangeInitialData()
{
  // F: send, receive, S: receive, send
  if (doesFirstStep()) {
    if (sendsInitializedData()) {
      sendData(getM2N());
    }
    if (receivesInitializedData()) {
      receiveData(getM2N());
    }
  } else { // second participant
    if (receivesInitializedData()) {
      receiveData(getM2N());
      // second participant has to save values for extrapolation
      if (isImplicitCouplingScheme() && getExtrapolationOrder() > 0) {
        doExtrapolationOn(getReceiveData());
      }
    }
    if (sendsInitializedData()) {
      if (isImplicitCouplingScheme() && getExtrapolationOrder() > 0) {
        doExtrapolationOn(getSendData());
      }
      sendData(getM2N());
    }
  }
}

std::pair<bool, bool> ParallelCouplingScheme::doAdvance()
{
  bool convergence, convergenceCoarseOptimization, doOnlySolverEvaluation; // @todo having the bools for convergence measurement declared for explicit and implicit coupling is not nice

  // initialize advance
  if (isImplicitCouplingScheme()) {
    convergenceCoarseOptimization = true;
    doOnlySolverEvaluation        = false;
  }

  // pre-acceleration communication
  if (doesFirstStep()) { //first participant
    PRECICE_DEBUG("Sending data...");
    sendData(getM2N());
  } else { //second participant
    PRECICE_DEBUG("Receiving data...");
    if (isExplicitCouplingScheme()) {
      receiveAndSetTimeWindowSize();
    }
    receiveData(getM2N());
  }

  // acceleration (only second participant)
  if (isImplicitCouplingScheme() && not doesFirstStep()) {
    ValuesMap designSpecifications; // TODO make this better?
    doAcceleration(designSpecifications, convergence, convergenceCoarseOptimization, doOnlySolverEvaluation);
  }

  // post-acceleration communication
  if (doesFirstStep()) { //first participant
    PRECICE_DEBUG("Receiving data...");
    if(isImplicitCouplingScheme()) {
      convergence = checkConvergence();
    } else {
      PRECICE_ASSERT(isExplicitCouplingScheme());
      receiveAndSetTimeWindowSize();
    }
    receiveData(getM2N());
  } else { //second participant
    PRECICE_DEBUG("Sending data...");
    if (isImplicitCouplingScheme()) {
      getM2N()->send(convergence);
      getM2N()->send(getIsCoarseModelOptimizationActive());
    }
    sendData(getM2N());
  }

  return std::pair<bool, bool>(convergence, convergenceCoarseOptimization);
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
