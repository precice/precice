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
                         secondParticipant, localParticipant, m2n, maxIterations, cplMode, dtMethod){}

void SerialCouplingScheme::initializeImplicit()
{
  PRECICE_CHECK(not getSendData().empty(), "No send data configured! Use explicit scheme for one-way coupling.");
  if (not doesFirstStep()) {
    if (not getConvergenceMeasures().empty()) {
      setupConvergenceMeasures();       // needs _couplingData configured
      setupDataMatrices(getSendData()); // Reserve memory and initialize data with zero
    }
    if (getAcceleration().get() != nullptr) {
      getAcceleration()->initialize(getSendData()); // Reserve memory, initialize
    }
  } else if (getAcceleration().get() != nullptr && not getAcceleration()->getDataIDs().empty()) {
    int dataID = *(getAcceleration()->getDataIDs().begin());
    PRECICE_CHECK(getSendData(dataID) == nullptr,
                  "In case of serial coupling, acceleration can be defined for "
                      << "data of second participant only!");
  }
}

void SerialCouplingScheme::initializeImplementation()
{
  for (DataMap::value_type &pair : getSendData()) {
    if (pair.second->initialize) {
      PRECICE_CHECK(not doesFirstStep(), "Only second participant can initialize data!");
      break;
    }
  }

  initializeSendingParticipants(getSendData());

  for (DataMap::value_type &pair : getReceiveData()) {
    if (pair.second->initialize) {
      PRECICE_CHECK(doesFirstStep(), "Only first participant can receive initial data!");
    }
  }

  initializeReceivingParticipants(getReceiveData());

  // If the second participant initializes data, the first receive for the
  // second participant is done in initializeData() instead of initialize().
  if (not doesFirstStep() && not hasToSendInitData() && isCouplingOngoing()) {
    PRECICE_DEBUG("Receiving data");
    receiveAndSetTimeWindowSize();
    receiveData(getM2N());
  }
}

void SerialCouplingScheme::exchangeInitialData()
{
  if (hasToReceiveInitData() && isCouplingOngoing()) {
    PRECICE_ASSERT(doesFirstStep());
    PRECICE_DEBUG("Receiving data");
    receiveData(getM2N());
  }

  if (hasToSendInitData() && isCouplingOngoing()) {
    PRECICE_ASSERT(not doesFirstStep());
    for (DataMap::value_type &pair : getSendData()) {
      if (pair.second->oldValues.cols() == 0)
        break;
      pair.second->oldValues.col(0) = *pair.second->values;
      // For extrapolation, treat the initial value as old time window value
      utils::shiftSetFirst(pair.second->oldValues, *pair.second->values);
    }

    // The second participant sends the initialized data to the first participant
    // here, which receives the data on call of initialize().
    sendData(getM2N());
    receiveAndSetTimeWindowSize();
    // This receive replaces the receive in initialize().
    receiveData(getM2N());
  }
}

void SerialCouplingScheme::explicitAdvance()
{
  PRECICE_DEBUG("Sending data...");
  sendTimeWindowSize();
  sendData(getM2N());

  if (isCouplingOngoing() || doesFirstStep()) {
    PRECICE_DEBUG("Receiving data...");
    receiveAndSetTimeWindowSize();
    receiveData(getM2N());
  }
}

std::pair<bool, bool> SerialCouplingScheme::implicitAdvance()
{
  bool convergence                   = true;
  bool convergenceCoarseOptimization = true;
  bool isCoarseModelOptimizationActive = false;
  bool doOnlySolverEvaluation        = false;

  PRECICE_DEBUG("Computed full length of iteration");
  if (doesFirstStep()) {
    sendTimeWindowSize();
    sendData(getM2N());
    getM2N()->receive(convergence);
    getM2N()->receive(isCoarseModelOptimizationActive);
    if(isCoarseModelOptimizationActive){
      activateCoarseModelOptimization();
    } else {
      deactivateCoarseModelOptimization();
    }
    if (convergence) {
      timeWindowCompleted();
    }
    receiveData(getM2N());
  } else {

    // get the current design specifications from the acceleration (for convergence measure)
    std::map<int, Eigen::VectorXd> designSpecifications;
    if (getAcceleration().get() != nullptr) {
      designSpecifications = getAcceleration()->getDesignSpecification(getSendData());
    }
    // measure convergence of coupling iteration
    // measure convergence for coarse model optimization
    if (getIsCoarseModelOptimizationActive()) {
      PRECICE_DEBUG("measure convergence of coarse model optimization.");
      // in case of multilevel acceleration only: measure the convergence of the coarse model optimization
      convergenceCoarseOptimization = measureConvergenceCoarseModelOptimization(designSpecifications);
      // Stop, when maximal iteration count (given in config) is reached
      if (maxIterationsReached())
        convergenceCoarseOptimization = true;

      convergence = false;
      // in case of multilevel PP only: if coarse model optimization converged
      // steering the requests for evaluation of coarse and fine model, respectively
      if (convergenceCoarseOptimization) {
        deactivateCoarseModelOptimization();
        doOnlySolverEvaluation           = true;
      } else {
        activateCoarseModelOptimization();
      }
    }
    // measure convergence of coupling iteration
    else {
      PRECICE_DEBUG("measure convergence.");
      doOnlySolverEvaluation = false;

      // measure convergence of the coupling iteration,
      convergence = measureConvergence(designSpecifications);
      // Stop, when maximal iteration count (given in config) is reached
      if (maxIterationsReached())
        convergence = true;
    }

    // for multi-level case, i.e., manifold mapping: after convergence of coarse problem
    // we only want to evaluate the fine model for the new input, no acceleration etc..
    if (not doOnlySolverEvaluation) {
      // coupling iteration converged for current time window. Advance in time.
      if (convergence) {
        if (getAcceleration().get() != nullptr) {
          setDeletedColumnsPPFiltering(getAcceleration()->getDeletedColumns());
          getAcceleration()->iterationsConverged(getSendData());
        }
        newConvergenceMeasurements();
        timeWindowCompleted();

        // no convergence achieved for the coupling iteration within the current time step
      } else if (getAcceleration().get() != nullptr) {
        getAcceleration()->performAcceleration(getSendData());
      }

      // extrapolate new input data for the solver evaluation in time.
      if (convergence && (getExtrapolationOrder() > 0)) {
        extrapolateData(getSendData()); // Also stores data
      } else {                          // Store data for conv. measurement, acceleration, or extrapolation
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

      /*
      /// @todo: (Edit: Done in the solver now) need to copy coarse old values to fine old values, as first solver always sends zeros to the second solver (as pressure vals)
      //       in the serial scheme, only the sendData is registered in MM PP, we also need to register the pressure values, i.e.
      //       old fine pressure vals = old coarse pressure vals TODO: find better solution,
      //auto fineIDs = getAcceleration()->getDataIDs();
      //for(auto id: fineIDs){
      //  std::cout<<"id: "<<id<<", fineIds.size(): "<<fineIDs.size()<<'\n';
      //  getReceiveData(id)->oldValues.column(0) = getReceiveData(id+fineIDs.size())->oldValues.column(0);
      //}
       */

      // only fine model solver evaluation is done, no PP
    } else {

      // if the coarse model problem converged within the first iteration, i.e., no acceleration at all
      // we need to register the coarse initialized data again on the fine input data,
      // otherwise the fine input data would be zero in this case, neither anything has been computed so far for the fine
      // model nor the acceleration did any data registration
      // ATTENTION: assumes that coarse data is defined after fine data in same ordering.
      if (getIterationsCoarseOptimization() == 1 && getAcceleration().get() != nullptr) {
        auto fineIDs = getAcceleration()->getDataIDs();
        for (auto &fineID : fineIDs) {
          (*getSendData(fineID)->values) = getSendData(fineID + fineIDs.size() + 1)->oldValues.col(0);
        }
      }
    }

    getM2N()->send(convergence);

    getM2N()->send(getIsCoarseModelOptimizationActive());

    sendData(getM2N());

    // the second participant does not want new data in the last iteration of the last time window
    if (isCouplingOngoing() || not convergence) {
      receiveAndSetTimeWindowSize();
      receiveData(getM2N());
    }
  }

  return std::pair<bool, bool>(convergence, convergenceCoarseOptimization);
}

} // namespace cplscheme
} // namespace precice
