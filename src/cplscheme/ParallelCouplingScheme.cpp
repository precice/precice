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
                         secondParticipant, localParticipant, m2n, maxIterations, cplMode, dtMethod){}

void ParallelCouplingScheme::initializeImplicit()
{
  PRECICE_CHECK(not getSendData().empty(), "No send data configured. Use explicit scheme for one-way coupling.");
  if (not doesFirstStep()) {         // second participant
    setupConvergenceMeasures();      // needs _couplingData configured
    mergeData();                     // merge send and receive data for all pp calls
    setupDataMatrices(getAcceleratedData()); // Reserve memory and initialize data with zero
    if (getAcceleration().get() != nullptr) {
      getAcceleration()->initialize(getAcceleratedData()); // Reserve memory, initialize
    }
  }
}

void ParallelCouplingScheme::initializeImplementation()
{
  if(anyDataRequiresInitialization(getSendData())) {
    hasToSendInitializedData();
  }

  if(anyDataRequiresInitialization(getReceiveData())) {
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
  }
  else { // second participant
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

void ParallelCouplingScheme::explicitAdvance()
{
  if (doesFirstStep()) {
    PRECICE_DEBUG("Sending data...");
    sendTimeWindowSize();
    sendData(getM2N());

    PRECICE_DEBUG("Receiving data...");
    receiveAndSetTimeWindowSize();
    receiveData(getM2N());
  } else { //second participant
    PRECICE_DEBUG("Receiving data...");
    receiveAndSetTimeWindowSize();
    receiveData(getM2N());

    PRECICE_DEBUG("Sending data...");
    sendTimeWindowSize();
    sendData(getM2N());
  }
}

std::pair<bool, bool> ParallelCouplingScheme::implicitAdvance()
{
  PRECICE_DEBUG("Computed full length of iteration");
  bool convergence, isCoarseModelOptimizationActive;
  bool convergenceCoarseOptimization   = true;
  bool doOnlySolverEvaluation          = false;
  if (doesFirstStep()) { //First participant
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
  } else { // second participant
    receiveData(getM2N());

    // get the current design specifications from the acceleration (for convergence measure)
    std::map<int, Eigen::VectorXd> designSpecifications;
    if (getAcceleration().get() != nullptr) {
      designSpecifications = getAcceleration()->getDesignSpecification(getAcceleratedData());
    }

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
      if (convergence) {
        if (getAcceleration().get() != nullptr) {
          setDeletedColumnsPPFiltering(getAcceleration()->getDeletedColumns());
          getAcceleration()->iterationsConverged(getAcceleratedData());
        }
        newConvergenceMeasurements();
        timeWindowCompleted();
      } else if (getAcceleration().get() != nullptr) {
        getAcceleration()->performAcceleration(getAcceleratedData());
      }

      // extrapolate new input data for the solver evaluation in time.
      if (convergence && (getExtrapolationOrder() > 0)) {
        extrapolateData(getAcceleratedData()); // Also stores data
      } else {                         // Store data for conv. measurement, acceleration, or extrapolation
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
    } else {

      // if the coarse model problem converged within the first iteration, i.e., no acceleration at all
      // we need to register the coarse initialized data again on the fine input data,
      // otherwise the fine input data would be zero in this case, neither anything has been computed so far for the fine
      // model nor the acceleration did any data registration
      // ATTENTION: assumes that coarse data is defined after fine data in same ordering.
      if (getIterationsCoarseOptimization() == 1 && getAcceleration().get() != nullptr) {
        auto  fineIDs         = getAcceleration()->getDataIDs();
        auto &acceleratedData = getAcceleratedData();
        for (auto &fineID : fineIDs) {
          *acceleratedData.at(fineID)->values = acceleratedData.at(fineID + fineIDs.size())->oldValues.col(0);
        }
      }
    }

    getM2N()->send(convergence);
    getM2N()->send(getIsCoarseModelOptimizationActive());

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
