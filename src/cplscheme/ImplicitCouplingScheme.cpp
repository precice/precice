// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "ImplicitCouplingScheme.hpp"
#include "impl/PostProcessing.hpp"
#include "impl/ConvergenceMeasure.hpp"
#include "Constants.hpp"
#include "mesh/SharedPointer.hpp"
#include "com/Communication.hpp"
#include "com/SharedPointer.hpp"
#include "io/TXTWriter.hpp"
#include "io/TXTReader.hpp"
#include "tarch/plotter/globaldata/TXTTableWriter.h"
#include <limits>

namespace precice {
namespace cplscheme {

tarch::logging::Log ImplicitCouplingScheme::
    _log("precice::cplscheme::ImplicitCouplingScheme" );

ImplicitCouplingScheme:: ImplicitCouplingScheme
(
  double                maxTime,
  int                   maxTimesteps,
  double                timestepLength,
  int                   validDigits,
  const std::string&    firstParticipant,
  const std::string&    secondParticipant,
  const std::string&    localParticipant,
  com::PtrCommunication communication,
  int                   maxIterations,
  constants::TimesteppingMethod dtMethod )
:
  BaseCouplingScheme(maxTime, maxTimesteps, timestepLength, validDigits),
  _firstParticipant(firstParticipant),
  _secondParticipant(secondParticipant),
  _doesFirstStep(false),
  _communication(communication),
  _iterationsWriter("iterations-" + localParticipant + ".txt"),
  //_residualWriterL1("residualL1-" + localParticipant + ".txt"),
  //_residualWriterL2("residualL2-" + localParticipant + ".txt"),
  //_amplificationWriter("amplification-" + localParticipant + ".txt"),
  _convergenceMeasures(),
  _postProcessing(),
  _extrapolationOrder(0),
  _maxIterations(maxIterations),
  _iterationToPlot(0),
  _timestepToPlot(0),
  _timeToPlot(0.0),
  _iterations(0),
  _totalIterations(0),
  _participantSetsDt(false),
  _participantReceivesDt(false)
{
  preciceCheck(_firstParticipant != _secondParticipant,
               "ImplicitCouplingScheme()", "First participant and "
               << "second participant must have different names!");
  if (dtMethod == constants::FIXED_DT){
    preciceCheck(not tarch::la::equals(timestepLength, UNDEFINED_TIMESTEP_LENGTH),
        "ImplicitCouplingScheme()", "Timestep length value has to be given "
        << "when the fixed timestep length method is chosen for an implicit "
        << "coupling scheme!");
  }
  if (localParticipant == _firstParticipant){
    _doesFirstStep = true;
    if (dtMethod == constants::FIRST_PARTICIPANT_SETS_DT){
      _participantSetsDt = true;
      setTimestepLength(UNDEFINED_TIMESTEP_LENGTH);
    }
  }
  else if (localParticipant == _secondParticipant){
    if (dtMethod == constants::FIRST_PARTICIPANT_SETS_DT){
      _participantReceivesDt = true;
    }
  }
  else {
    preciceError("initialize()", "Name of local participant \""
                 << localParticipant << "\" does not match any "
                 << "participant specified for the coupling scheme!");
  }
  preciceCheck((maxIterations > 0) || (maxIterations == -1),
               "ImplicitCouplingState()",
               "Maximal iteration limit has to be larger than zero!");
  assertion(_communication.use_count() > 0);
}

ImplicitCouplingScheme:: ~ImplicitCouplingScheme()
{}

void ImplicitCouplingScheme:: setExtrapolationOrder
(
  int order )
{
  preciceCheck((order == 0) || (order == 1) || (order == 2),
               "setExtrapolationOrder()", "Extrapolation order has to be "
               << " 0, 1, or 2!");
  _extrapolationOrder = order;
}

void ImplicitCouplingScheme:: addConvergenceMeasure
(
  int                         dataID,
  bool                        suffices,
  impl::PtrConvergenceMeasure measure )
{
  ConvergenceMeasure convMeasure;
  convMeasure.dataID = dataID;
  convMeasure.data = NULL;
  convMeasure.suffices = suffices;
  convMeasure.measure = measure;
  _convergenceMeasures.push_back(convMeasure);
}

void ImplicitCouplingScheme:: setIterationPostProcessing
(
  impl::PtrPostProcessing postProcessing )
{
  assertion(postProcessing.get() != NULL);
  _postProcessing = postProcessing;
}

void ImplicitCouplingScheme:: initialize
(
  double startTime,
  int    startTimestep )
{
  preciceTrace2("initialize()", startTime, startTimestep);
  assertion(not isInitialized());
  assertion1(tarch::la::greaterEquals(startTime, 0.0), startTime);
  assertion1(startTimestep >= 0, startTimestep);
  assertion(_communication->isConnected());
  preciceCheck(not getSendData().empty(), "initialize()",
               "No send data configured! Use explicit scheme for one-way coupling.");
  setTime(startTime);
  setTimesteps(startTimestep);
  if (not _doesFirstStep){
    setupConvergenceMeasures(); // needs _couplingData configured
    setupDataMatrices(); // Reserve memory and initialize data with zero
    if (_postProcessing.get() != NULL){
      _postProcessing->initialize(getSendData()); // Reserve memory, initialize
    }
  }
  else if (_postProcessing.get() != NULL){
    int dataID = _postProcessing->getDataID();
    preciceCheck(getSendData(dataID) == NULL, "initialize()",
                 "A post-processing can be defined for data of second "
                 << "participant only!");
  }

  requireAction(constants::actionWriteIterationCheckpoint());

  // Determine data initialization
  bool doesReceiveData = not _doesFirstStep;

  // If the second participant initializes data, the first receive for the
  // second participant is done in initializeData() instead of initialize().
  foreach (DataMap::value_type & pair, getSendData()){
    if (pair.second.initialize){
      preciceCheck(not _doesFirstStep, "initialize()",
                   "Only second participant can initialize data!");
      requireAction(constants::actionWriteInitialData());
      preciceDebug("Initialized data to be written");
      doesReceiveData = false;
      break;
    }
  }
  // If the second participant initializes data, the first receive for the first
  // participant is done in initialize() instead of andvance().
  foreach (DataMap::value_type & pair, getReceiveData()){
    if (pair.second.initialize){
      preciceCheck(_doesFirstStep, "initialize()",
                   "Only first participant can receive initial data!");
      preciceDebug("Initialized data to be received");
      doesReceiveData = true;
    }
  }

  if (doesReceiveData && isCouplingOngoing()){
    preciceDebug("Receiving data");
    _communication->startReceivePackage(0);
    if (_participantReceivesDt){
      double dt = UNDEFINED_TIMESTEP_LENGTH;
      _communication->receive(dt, 0);
      preciceDebug("received timestep length of " << dt);
      assertion(not tarch::la::equals(dt, UNDEFINED_TIMESTEP_LENGTH));
      setTimestepLength(dt);
      //setMaxLengthNextTimestep(dt);
    }
    receiveData(_communication);
    _communication->finishReceivePackage();
    setHasDataBeenExchanged(true);
  }
  initializeTXTWriters();
  setIsInitialized(true);
}

void ImplicitCouplingScheme:: initializeData()
{
  preciceTrace("initializeData()");
  preciceCheck(isInitialized(), "initializeData()",
               "initializeData() can be called after initialize() only!");
  preciceCheck(isActionRequired(constants::actionWriteInitialData()),
               "initializeData()", "Not required data initialization!");
  assertion(not _doesFirstStep);
  foreach (DataMap::value_type & pair, getSendData()){
    utils::DynVector& oldValues = pair.second.oldValues.column(0);
    oldValues = *pair.second.values;

    // For extrapolation, treat the initial value as old timestep value
    pair.second.oldValues.shiftSetFirst(*pair.second.values);

    // The second participant sends the initialized data to the first particpant
    // here, which receives the data on call of initialize().
    sendData(_communication);
    _communication->startReceivePackage(0);
    // This receive replaces the receive in initialize().
    receiveData(_communication);
    _communication->finishReceivePackage();
    setHasDataBeenExchanged(true);
  }
  performedAction(constants::actionWriteInitialData());
}

//void ImplicitCouplingScheme:: addComputedTime
//(
//  double timeToAdd )
//{
//  preciceTrace2("addComputedTime()", timeToAdd, getTime());
//  preciceCheck(isCouplingOngoing(), "addComputedTime()",
//               "Invalid call of addComputedTime() after simulation end!");
//
//  // Check validness
//  double eps = std::pow(10.0, -1 * getValidDigits());
//  bool greaterThanZero = tarch::la::greater(timeToAdd, 0.0, eps);
//  preciceCheck(greaterThanZero, "addComputedTime()", "The computed timestep length "
//               << "exceeds the maximum timestep limit for this time step!");
//
//  setComputedTimestepPart(getComputedTimestepPart() + timeToAdd);
//  setTime(getTime() + timeToAdd);
//  //setSubIteration(_iterations + 1);
//  //_totalIterations++;
//}

void ImplicitCouplingScheme:: advance()
{
  preciceTrace2("advance()", getTimesteps(), getTime());
  checkCompletenessRequiredActions();
  setHasDataBeenExchanged(false);
  setIsCouplingTimestepComplete(false);
  double eps = std::pow(10.0, -1 * getValidDigits());
  bool convergence = false;
  if (tarch::la::equals(getThisTimestepRemainder(), 0.0, eps)){
    preciceDebug("Computed full length of iteration");
    if (_doesFirstStep){
      _communication->startSendPackage(0);
      if (_participantSetsDt){
        preciceDebug("sending timestep length of " << getComputedTimestepPart());
        _communication->send(getComputedTimestepPart(), 0);
      }
      sendData(_communication);
      _communication->finishSendPackage();
      _communication->startReceivePackage(0);
      _communication->receive(convergence, 0);
      if (convergence){
        timestepCompleted();
      }
      if (isCouplingOngoing()){
        receiveData(_communication);
      }
      _communication->finishReceivePackage();
    }
    else {
      convergence = measureConvergence();
      assertion2((_iterations <= _maxIterations) || (_maxIterations == -1),
                 _iterations, _maxIterations);
      // Stop, when maximal iteration count (given in config) is reached
      if (_iterations == _maxIterations-1){
        convergence = true;
      }
      if (convergence){
        if (_postProcessing.get() != NULL){
          _postProcessing->iterationsConverged(getSendData());
        }
        newConvergenceMeasurements();
        timestepCompleted();
      }
      else if (_postProcessing.get() != NULL){
        _postProcessing->performPostProcessing(getSendData());
      }
      _communication->startSendPackage(0);
      _communication->send(convergence, 0);
      if (isCouplingOngoing()){
        if (convergence && (_extrapolationOrder > 0)){
          extrapolateData(); // Also stores data
        }
        else { // Store data for conv. measurement, post-processing, or extrapolation
          foreach (DataMap::value_type& pair, getSendData()){
            if (pair.second.oldValues.size() > 0){
              pair.second.oldValues.column(0) = *pair.second.values;
            }
          }
          foreach (DataMap::value_type& pair, getReceiveData()){
            if (pair.second.oldValues.size() > 0){
              pair.second.oldValues.column(0) = *pair.second.values;
            }
          }
        }
        sendData(_communication);
        _communication->finishSendPackage();
        _communication->startReceivePackage(0);
        if (_participantReceivesDt){
          double dt = UNDEFINED_TIMESTEP_LENGTH;
          _communication->receive(dt, 0);
          assertion(not tarch::la::equals(dt, UNDEFINED_TIMESTEP_LENGTH));
          setTimestepLength(dt);
        }
        receiveData(_communication);
        _communication->finishReceivePackage();
      }
      else {
        _communication->finishSendPackage();
      }
    }

    if (not convergence){
      preciceDebug("No convergence achieved");
      requireAction(constants::actionReadIterationCheckpoint());
      _iterations++;
      _totalIterations++;
      // The computed timestep part equals the timestep length, since the
      // timestep remainder is zero. Subtract the timestep length do another
      // coupling iteration.
      assertion(tarch::la::greater(getComputedTimestepPart(), 0.0));
      setTime(getTime() - getComputedTimestepPart());
    }
    else {
      preciceDebug("Convergence achieved");
      _iterationsWriter.writeData("Timesteps", getTimesteps());
      _iterationsWriter.writeData("Total Iterations", _totalIterations);
      _iterationsWriter.writeData("Iterations", _iterations);
      int converged = _iterations < _maxIterations ? 1 : 0;
      _iterationsWriter.writeData("Convergence", converged);
      _iterations = 0;
    }
    setHasDataBeenExchanged(true);
    setComputedTimestepPart(0.0);
  }

  // When the iterations of one timestep are converged, the old time, timesteps,
  // and iteration should be plotted, and not the 0th of the new timestep. Thus,
  // the plot values are only updated when no convergence was achieved.
  if (not convergence){
    _timestepToPlot = getTimesteps();
    _timeToPlot = getTime();
    _iterationToPlot = _iterations;
  }
  else {
    _iterationToPlot++;
  }
}

void ImplicitCouplingScheme:: timestepCompleted()
{
  preciceTrace2("timestepCompleted()", getTimesteps(), getTime());
  preciceInfo("timestepCompleted()", "Timestep completed");
  setIsCouplingTimestepComplete(true);
  setTimesteps(getTimesteps() + 1 );
  //setTime(getTimesteps() * getTimestepLength() ); // Removes numerical errors
  if (isCouplingOngoing()){
    preciceDebug("Setting require create checkpoint");
    requireAction(constants::actionWriteIterationCheckpoint());
  }
}

void ImplicitCouplingScheme:: finalize()
{
   preciceTrace("finalize()");
   checkCompletenessRequiredActions();
   preciceCheck(isInitialized(), "finalize()",
                "Called finalize() before initialize()!");
   preciceCheck(not isCouplingOngoing(), "finalize()",
                "Called finalize() while isCouplingOngoing() returns true!");
}

void ImplicitCouplingScheme:: initializeTXTWriters()
{
  _iterationsWriter.addData("Timesteps", io::TXTTableWriter::INT );
  _iterationsWriter.addData("Total Iterations", io::TXTTableWriter::INT );
  _iterationsWriter.addData("Iterations", io::TXTTableWriter::INT );
  _iterationsWriter.addData("Convergence", io::TXTTableWriter::INT );
//  _residualWriterL1.addData("Iterations", io::TXTTableWriter::INT );
//  _residualWriterL2.addData ("Iterations", io::TXTTableWriter::INT );
//  _amplificationWriter.addData("Iterations", io::TXTTableWriter::INT );
//  _residualWriterL1.addData("L1-Residual", io::TXTTableWriter::DOUBLE );
//  _residualWriterL2.addData("L2-Residual", io::TXTTableWriter::DOUBLE );
//  _amplificationWriter.addData("Amplification" ,io::TXTTableWriter::DOUBLE );
//  int entries = 0;
//  if(getSendData().size() > 0 ){
//    entries = (int)getSendData().begin()->second.values->size();
//  }
//  int levels = 1;
//  int treatedEntries = 2;
//  int entriesCurrentLevel = 1;
//  while (treatedEntries < entries){
//    treatedEntries += entriesCurrentLevel;
//    levels ++;
//    entriesCurrentLevel *= 2;
//  }
//  if (treatedEntries == entries){
//    for (int i=0; i < levels; i++ ){
//      _residualWriterL1.addData("L1-Residual-level-"+i, io::TXTTableWriter::DOUBLE);
//      _residualWriterL2.addData("L2-Residual-level-"+i, io::TXTTableWriter::DOUBLE);
//      _amplificationWriter.addData("Amplification-level-"+i,io::TXTTableWriter::DOUBLE);
//    }
//  }
}

void ImplicitCouplingScheme:: setupDataMatrices()
{
  preciceTrace("setupDataMatrices()");
  // Reserve storage for convergence measurement of send and receive data values
  foreach (ConvergenceMeasure& convMeasure, _convergenceMeasures){
    assertion(convMeasure.data != NULL);
    if (convMeasure.data->oldValues.cols() < 1){
      convMeasure.data->oldValues.append(CouplingData::DataMatrix(
          convMeasure.data->values->size(), 1, 0.0));
    }
  }
  // Reserve storage for extrapolation of send data values
  if (_extrapolationOrder > 0){
    foreach (DataMap::value_type& pair, getSendData()){
      int cols = pair.second.oldValues.cols();
      assertion1(cols <= 1, cols);
      pair.second.oldValues.append(CouplingData::DataMatrix(
          pair.second.values->size(), _extrapolationOrder + 1 - cols, 0.0));
    }
  }
}

void ImplicitCouplingScheme:: setupConvergenceMeasures()
{
  preciceTrace("setupConvergenceMeasures()");
  assertion(not _doesFirstStep);
  preciceCheck(not _convergenceMeasures.empty(), "setupConvergenceMeasures()",
      "At least one convergence measure has to be defined for "
      << "an implicit coupling scheme!");
  foreach (ConvergenceMeasure& convMeasure, _convergenceMeasures){
    int dataID = convMeasure.dataID;
    if ((getSendData(dataID) != NULL)){
      convMeasure.data = getSendData(dataID);
    }
    else {
      convMeasure.data = getReceiveData(dataID);
      assertion(convMeasure.data != NULL);
    }
  }
}

bool ImplicitCouplingScheme:: measureConvergence()
{
  preciceTrace("measureLocalConvergence()");
  using boost::get;
  bool allConverged = true;
  bool oneSuffices = false;
  assertion(_convergenceMeasures.size() > 0);
  foreach(ConvergenceMeasure& convMeasure, _convergenceMeasures){
    assertion(convMeasure.data != NULL);
    assertion(convMeasure.measure.get() != NULL);
    utils::DynVector& oldValues = convMeasure.data->oldValues.column(0);
    convMeasure.measure->measure(oldValues, *convMeasure.data->values);
    if (not convMeasure.measure->isConvergence()){
      //preciceDebug("Local convergence = false");
      allConverged = false;
    }
    else if (convMeasure.suffices == true){
      oneSuffices = true;
    }
    preciceInfo("measureConvergence()", convMeasure.measure->printState());
  }
  if (allConverged){
    preciceInfo("measureConvergence()", "All converged");
  }
  else if (oneSuffices){
    preciceInfo("measureConvergence()", "Sufficient measure converged");
  }
  return allConverged || oneSuffices;
}

void ImplicitCouplingScheme:: extrapolateData()
{
   preciceTrace("extrapolateData()");
   bool startWithFirstOrder = (getTimesteps() == 1) && (_extrapolationOrder == 2);
   if((_extrapolationOrder == 1) || startWithFirstOrder ){
      preciceInfo("extrapolateData()", "Performing first order extrapolation" );
      foreach(DataMap::value_type & pair, getSendData() ){
         assertion(pair.second.oldValues.cols() > 1 );
         utils::DynVector & values = *pair.second.values;
         pair.second.oldValues.column(0) = values;    // = x^t
         values *= 2.0;                                  // = 2 * x^t
         values -= pair.second.oldValues.column(1);   // = 2*x^t - x^(t-1)
         pair.second.oldValues.shiftSetFirst(values );
      }
   }
   else if(_extrapolationOrder == 2 ){
      preciceInfo("extrapolateData()", "Performing second order extrapolation" );
      foreach(DataMap::value_type & pair, getSendData() ) {
         assertion(pair.second.oldValues.cols() > 2 );
         utils::DynVector & values = *pair.second.values;
         pair.second.oldValues.column(0) = values;        // = x^t                                     // = 2.5 x^t
//         utils::DynVector & valuesOld1 = pair.second.oldValues.getColumn(1);
//         utils::DynVector & valuesOld2 = pair.second.oldValues.getColumn(2);
//         for(int i=0; i < values.size(); i++ ) {
//            values[i] -= valuesOld1[i] * 3.0; // =
//            values[i] += valuesOld2[i] * 3.0; // =
//         }
         values *= 2.5;                                      // = 2.5 x^t
         utils::DynVector & valuesOld1 = pair.second.oldValues.column(1);
         utils::DynVector & valuesOld2 = pair.second.oldValues.column(2);
         for(int i=0; i < values.size(); i++ ){
            values[i] -= valuesOld1[i] * 2.0; // = 2.5x^t - 2x^(t-1)
            values[i] += valuesOld2[i] * 0.5; // = 2.5x^t - 2x^(t-1) + 0.5x^(t-2)
         }
         pair.second.oldValues.shiftSetFirst(values );
         //preciceDebug("extrapolateData()", "extrapolated data to \""
         //               << *pair.second.values );
      }
   }
   else {
      preciceError("extrapolateData()", "Called extrapolation with order != 1,2!" );
   }
}

void ImplicitCouplingScheme:: newConvergenceMeasurements()
{
   preciceTrace("newConvergenceMeasurements()");
   foreach (ConvergenceMeasure& convMeasure, _convergenceMeasures){
     assertion(convMeasure.measure.get() != NULL);
     convMeasure.measure->newMeasurementSeries();
   }
}

std::vector<std::string> ImplicitCouplingScheme:: getCouplingPartners() const
{
  std::vector<std::string> partnerNames;

  if(_doesFirstStep){
    partnerNames.push_back(_firstParticipant);
  }
  else {
    partnerNames.push_back(_secondParticipant);
  }
  return partnerNames;
}

void ImplicitCouplingScheme:: sendState
(
  com::PtrCommunication communication,
  int                   rankReceiver )
{
  preciceTrace1("sendState()", rankReceiver);
  communication->startSendPackage(rankReceiver );
  BaseCouplingScheme::sendState(communication, rankReceiver );
  communication->send(_maxIterations, rankReceiver );
  communication->send(_iterations, rankReceiver );
  communication->send(_totalIterations, rankReceiver );
  communication->finishSendPackage();
}

void ImplicitCouplingScheme:: receiveState
(
  com::PtrCommunication communication,
  int                   rankSender )
{
  preciceTrace1("receiveState()", rankSender);
  communication->startReceivePackage(rankSender);
  BaseCouplingScheme::receiveState(communication, rankSender);
  communication->receive(_maxIterations, rankSender);
  int subIteration = -1;
  communication->receive(subIteration, rankSender);
  _iterations = subIteration;
  communication->receive(_totalIterations, rankSender);
  communication->finishReceivePackage();
}

std::string ImplicitCouplingScheme:: printCouplingState() const
{
  std::ostringstream os;
  os << " it " << _iterationToPlot; //_iterations;
  if(_maxIterations != -1 ){
    os << " of " << _maxIterations;
  }
  os << " | " << printBasicState(_timestepToPlot, _timeToPlot) << std::endl << printActionsState();
  return os.str();
}

void ImplicitCouplingScheme:: exportState
(
  const std::string& filenamePrefix ) const
{
  if (not _doesFirstStep){
    io::TXTWriter writer(filenamePrefix + "_cplscheme.txt");
    foreach (const BaseCouplingScheme::DataMap::value_type& dataMap, getSendData()){
      writer.write(dataMap.second.oldValues);
    }
    foreach (const BaseCouplingScheme::DataMap::value_type& dataMap, getReceiveData()){
      writer.write(dataMap.second.oldValues);
    }
    if (_postProcessing.get() != NULL){
      _postProcessing->exportState(writer);
    }
  }
}

void ImplicitCouplingScheme:: importState
(
  const std::string& filenamePrefix )
{
  if (not _doesFirstStep){
    io::TXTReader reader(filenamePrefix + "_cplscheme.txt");
    foreach (BaseCouplingScheme::DataMap::value_type& dataMap, getSendData()){
      reader.read(dataMap.second.oldValues);
    }
    foreach (BaseCouplingScheme::DataMap::value_type& dataMap, getReceiveData()){
      reader.read(dataMap.second.oldValues);
    }
    if (_postProcessing.get() != NULL){
      _postProcessing->importState(reader);
    }
  }
}

//void ImplicitCouplingScheme:: writeResidual
//(
//  const utils::DynVector& values,
//  const utils::DynVector& oldValues )
//{
//  using namespace tarch::la;
//  size_t entries = (size_t)values.size();
//
//  // Compute current residuals
//  utils::DynVector residual(values );
//  residual -= oldValues;
//  utils::DynVector oldValuesTemp(oldValues );
//
//  size_t levels = 1;
//  size_t treatedEntries = 2;
//  size_t entriesCurrentLevel = 1;
//  while(treatedEntries < entries ){
//    treatedEntries += entriesCurrentLevel;
//    levels ++;
//    entriesCurrentLevel *= 2;
//  }
//
//  double residualNormL1 = norm1(residual );
//  double residualNormL2 = norm2(residual );
//  double amplification = norm2(values) / norm2(oldValues);
//  utils::DynVector hierarchicalResidualNormsL1;
//  utils::DynVector hierarchicalResidualNormsL2;
//  utils::DynVector hierarchicalAmplification;
//
//  _residualWriterL1.writeData("L1-Residual", residualNormL1 );
//  _residualWriterL2.writeData("L2-Residual", residualNormL2 );
//  _amplificationWriter.writeData("Amplification", amplification );
//
//  if(treatedEntries == entries ){ // Hierarchizable
//    treatedEntries = 2;
//    hierarchicalResidualNormsL1.append(utils::DynVector(levels, 0.0) );
//    hierarchicalResidualNormsL2.append(utils::DynVector(levels, 0.0) );
//    hierarchicalAmplification.append(utils::DynVector(levels, 0.0) );
//    entriesCurrentLevel = std::pow(2.0, (int)(levels - 2));
//    for(size_t level=levels-1; level > 0; level-- ){
//      size_t stepsize = (entries - 1) / std::pow(2.0, (int)(level-1));
//      assertion(stepsize % 2 == 0 );
//      size_t index = stepsize / 2;
//
//      double amplificationNom = 0.0;
//      double amplificationDenom = 0.0;
//      for(size_t i=0; i < entriesCurrentLevel; i++ ){
//        residual[index] -=(residual[index - stepsize/2] +
//                             residual[index + stepsize/2] ) / 2.0;
//        oldValuesTemp[index] -=(oldValuesTemp[index - stepsize/2] +
//                                  oldValuesTemp[index + stepsize/2] ) / 2.0;
//        hierarchicalResidualNormsL1[level] += std::abs(residual[index] );
//        hierarchicalResidualNormsL2[level] += residual[index] * residual[index];
//        amplificationNom += std::pow(residual[index] + oldValuesTemp[index], 2 );
//        amplificationDenom += std::pow(oldValuesTemp[index], 2 );
//        index += stepsize;
//      }
//      hierarchicalResidualNormsL2[level] = std::sqrt(hierarchicalResidualNormsL2[level] );
//      hierarchicalAmplification[level] = std::sqrt(amplificationNom) / std::sqrt(amplificationDenom);
//      treatedEntries += entriesCurrentLevel;
//      entriesCurrentLevel /= 2;
//    }
//    hierarchicalResidualNormsL1[0] =
//        std::abs(residual[0]) + std::abs(residual[entries-1] );
//    hierarchicalResidualNormsL2[0] = std::sqrt(
//        residual[0] * residual[0] + residual[entries-1] * residual[entries-1] );
//    hierarchicalAmplification[0] =
//        std::sqrt(std::pow(values[0], 2) + std::pow(values[entries-1], 2)) /
//        std::sqrt(std::pow(oldValues[0], 2) + std::pow(oldValues[entries-1], 2));
//    for(size_t i=0; i < levels; i++ ){
//      _residualWriterL1.writeData("L1-Residual-level-" + i, residualNormL1 );
//      _residualWriterL2.writeData("L2-Residual-level-" + i, residualNormL2 );
//      _amplificationWriter.writeData("Amplification-level-" + i, amplification );
//    }
//  }
//}

}} // namespace precice, cplscheme
