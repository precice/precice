// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_CPLSCHEME_IMPLICITCOUPLINGSCHEME_HPP_
#define PRECICE_CPLSCHEME_IMPLICITCOUPLINGSCHEME_HPP_

#include "BaseCouplingScheme.hpp"
#include "SharedPointer.hpp"
#include "Constants.hpp"
#include "impl/SharedPointer.hpp"
#include "io/TXTTableWriter.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/PropertyContainer.hpp"
#include "mesh/Mesh.hpp"
#include "tarch/logging/Log.h"
#include "utils/Helpers.hpp"
#include "tarch/la/DynamicColumnMatrix.h"
#include "boost/tuple/tuple.hpp"

namespace precice {
  namespace cplscheme {
    namespace tests {
      class ImplicitCouplingSchemeTest;
    }
  }
}

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace cplscheme {

/**
 * TODO Abstract class, add to docu.
 *
 * @brief Coupling scheme with iterations per timestep to achieve strong solution.
 */
class ImplicitCouplingScheme : public BaseCouplingScheme
{
public:

  /**
   * @brief Constructor.
   *
   * TODO adapt params docu, e.g., first and second participant
   *
   * @param maxTime [IN] Simulation time limit, or UNDEFINED_TIME.
   * @param maxTimesteps [IN] Simulation timestep limit, or UNDEFINED_TIMESTEPS.
   * @param timestepLength [IN] Simulation timestep length.
   * @param firstParticipant [IN] Name of participant starting simulation.
   * @param secondParticipant [IN] Name of second participant in coupling.
   * @param localParticipant [IN] Name of participant using this coupling scheme.
   * @param communication [IN] Communication object for com. between participants.
   * @param maxIterations [IN] Maximal iterations per coupling timestep.
   * @param monitorIterations [IN] If true, a txt file monitoring iterations is
   *                          written.
   */
  ImplicitCouplingScheme (
    double                maxTime,
    int                   maxTimesteps,
    double                timestepLength,
    int                   validDigits,
    const std::string&    firstParticipant,
    const std::string&    secondParticipant,
    const std::string&    localParticipant,
    com::PtrCommunication communication,
    int                   maxIterations,
    constants::TimesteppingMethod dtMethod);

  /**
   * @brief Destructor.
   */
  virtual ~ImplicitCouplingScheme();

  /**
   * @brief Sets order of predictor of interface values for first participant.
   *
   * The first participant in the implicit coupling scheme has to take some
   * initial guess for the interface values computed by the second participant.
   * In order to improve this initial guess, an extrapolation from previous
   * timesteps can be performed.
   *
   * The standard predictor is of order zero, i.e., simply the converged values
   * of the last timestep are taken as initial guess for the coupling iterations.
   * Currently, an order 1 predictor is implement besides that.
   */
  void setExtrapolationOrder ( int order );

  /**
   * @brief Adds a measure to determine the convergence of coupling iterations.
   */
  void addConvergenceMeasure (
    int                         dataID,
    bool                        suffices,
    impl::PtrConvergenceMeasure measure );

  /**
   * @brief Set a coupling iteration post-processing technique.
   */
  void setIterationPostProcessing ( impl::PtrPostProcessing postProcessing );

  /**
   * @brief Initializes the coupling scheme.
   */
  virtual void initialize (
    double startTime,
    int    startTimestep ) =0;

  /**
   * @brief Initializes data with written values.
   *
   * Preconditions:
   * - initialize() has been called.
   * - advance() has NOT yet been called.
   */
  virtual void initializeData() =0;


  /**
   * @brief Advances within the coupling scheme (not necessarily in time).
   *
   * Preconditions:
   * - initialize() has been called.
   */
  virtual void advance() =0;

  /**
   * @brief Finalizes the coupling scheme.
   */
  virtual void finalize();

  /*
   * @brief returns list of all coupling partners
   */
  virtual std::vector<std::string> getCouplingPartners () const;

  virtual void sendState (
    com::PtrCommunication communication,
    int                   rankReceiver );

  virtual void receiveState (
    com::PtrCommunication communication,
    int                   rankSender );

  virtual std::string printCouplingState() const;

  virtual void exportState(const std::string& filenamePrefix) const;

  virtual void importState(const std::string& filenamePrefix);

protected:

  // @return True, if local participant is the one starting the explicit scheme.
  bool doesFirstStep(){
    return _doesFirstStep;
  }

  // @return Communication device to the other coupling participant.
  com::PtrCommunication getCommunication(){
    return _communication;
  }

  // TODO Refactor to getter/setter/other method or make private
  // @brief Responsible for monitoring iteration count over timesteps.
  io::TXTTableWriter _iterationsWriter;

  void setIterationToPlot(int iterationToPlot){
    _iterationToPlot = iterationToPlot;
  }

  void setTimestepToPlot(int timestepToPlot){
    _timestepToPlot = timestepToPlot;
  }

  void setTimeToPlot(double timeToPlot){
    _timeToPlot = timeToPlot;
  }

  // @return Post-processing method to speedup iteration convergence.
  impl::PtrPostProcessing getPostProcessing(){
    return _postProcessing;
  }

  void setHasToSendInitData(bool hasToSendInitData){
    _hasToSendInitData = hasToSendInitData;
  }

  void setHasToReceiveInitData(bool hasToReceiveInitData){
    _hasToReceiveInitData = hasToReceiveInitData;
  }

  bool hasToSendInitData(){
    return _hasToSendInitData;
  }

  bool hasToReceiveInitData(){
    return _hasToReceiveInitData;
  }

  bool participantReceivesDt(){
    return _participantReceivesDt;
  }

  bool participantSetsDt(){
    return _participantSetsDt;
  }

  void setIterations(int iterations){
    _iterations = iterations;
  }

  int getIterations(){
    return _iterations;
  }

  int getTotalIterations(){
    return _totalIterations;
  }

  void increaseIterations(){
    _iterations++;
  }

  void increaseTotalIterations(){
    _totalIterations++;
  }

  void increaseIterationToPlot(){
    _iterationToPlot++;
  }

  int getMaxIterations(){
    return _maxIterations;
  }

  int getExtrapolationOrder(){
    return _extrapolationOrder;
  }

  /**
   * @brief Initializes the txt writers for writing residuals, iterations, ...
   */
  void initializeTXTWriters();

  /**
   * @brief Sets up _dataStorage to store data values of last timestep.
   *
   * Every send data has an entry in _dataStorage. Every Entry is a vector
   * of data values with length according to the total number of values on all
   * meshes. The ordering of the data values corresponds to that in the meshes
   * and the ordering of the meshes to that in _couplingData.
   */
  void setupDataMatrices(DataMap& data);

  void setupConvergenceMeasures();

  void newConvergenceMeasurements();

  /**
   * @brief Updates internal state of coupling scheme for next timestep.
   */
  void timestepCompleted();

  /**
   * @brief Updates the convergence measurement of local send data.
   */
  bool measureConvergence();

  void extrapolateData(DataMap& data);

private:

  // @brief True, if local participant is the one starting the explicit scheme.
    bool _doesFirstStep;

  // @brief Communication device to the other coupling participant.
    com::PtrCommunication _communication;

  /**
   * @brief Holds relevant variables to perform a convergence measurement.
   */
  struct ConvergenceMeasure
  {
    int dataID;
    CouplingData* data;
    bool suffices;
    impl::PtrConvergenceMeasure measure;
  };

  typedef tarch::la::DynamicColumnMatrix<double> DataMatrix;

  typedef tarch::la::DynamicVector<double> DataVector;

  // @brief Logging device.
  static tarch::logging::Log _log;

  // @brief First participant name.
  std::string _firstParticipant;

  // @brief Second participant name.
  std::string _secondParticipant;

  // @brief Writes residuals to file.
//  io::TXTTableWriter _residualWriterL1;
//  io::TXT_communicationTableWriter _residualWriterL2;

  // @brief Writes value amplification to file.
//  io::TXTTableWriter _amplificationWriter;

  typedef boost::tuple<int,impl::PtrConvergenceMeasure> MeasureTuple;

  // @brief All convergence measures of coupling iterations.
  //
  // Before initialization, only dataID and measure variables are filled. Then,
  // the data is fetched from send and receive data assigned to the cpl scheme.
  std::vector<ConvergenceMeasure> _convergenceMeasures;

  // @brief Post-processing method to speedup iteration convergence.
  impl::PtrPostProcessing _postProcessing;

  // @brief Extrapolation order of coupling data for first iteration of every dt.
  int _extrapolationOrder;

  // @brief Limit of iterations during one timestep.
  int _maxIterations;

  // @brief Number of iteration in current timestep.
  int _iterationToPlot;
  int _timestepToPlot;
  double _timeToPlot;

  // @brief Number of iterations in current timestep.
  int _iterations;

  // @brief Number of total iterations performed.
  int _totalIterations;

  // @brief Determines, if the timestep length is set by the participant.
  bool _participantSetsDt;

  // @brief Determines, if the dt length is set received from the other participant
  bool _participantReceivesDt;

  // @brief to carry initData information from initialize to initData
  bool _hasToSendInitData;

  // @brief to carry initData information from initialize to initData
  bool _hasToReceiveInitData;


//  void writeResidual (
//    const utils::DynVector& values,
//    const utils::DynVector& oldValues );

  friend class tests::ImplicitCouplingSchemeTest;
};

}} // namespace precice, cplscheme

#endif /* PRECICE_CPLSCHEME_IMPLICITCOUPLINGSCHEME_HPP_ */
